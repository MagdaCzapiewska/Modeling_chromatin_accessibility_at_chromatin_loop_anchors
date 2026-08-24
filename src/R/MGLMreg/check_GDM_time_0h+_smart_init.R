library(data.table)

config_path <- file.path("config", "config.yml")
config <- yaml::yaml.load_file(config_path)

datadir <- config$paths$datadir
resultsdir <- config$paths$resultsdir

input_dir <- file.path(resultsdir, "MGLMreg", "GDM_time_0h+")
output_dir <- file.path(resultsdir, "MGLMreg")
output_file <- file.path(output_dir, "summary_GDM_time_0h+_smart_init.tsv")

loops_file <- file.path(datadir, "long_and_short_range_loops_D_mel.tsv")

loops_table <- fread(loops_file)
loops_of_interest <- loops_table$loop_id
loops_of_interest <- loops_of_interest[order(as.integer(sub("^L", "", loops_of_interest)))]

get_coef_safe <- function(mat, row, col) {
  if (!is.null(mat) && row %in% rownames(mat) && col %in% colnames(mat)) {
    return(mat[row, col])
  } else {
    return(NA_real_)
  }
}

results_list <- list()

for (loop_id in loops_of_interest) {
  
  cat("Loop:", loop_id, "\n")
  file <- file.path(input_dir, paste0("fit_reg_", loop_id, ".rds"))
  
  # Flaga czy udało się wczytać poprawny model
  is_model_valid <- FALSE
  if (!file.exists(file)) {
    next
  }
  
  if (file.exists(file)) {
    fit_reg <- readRDS(file)
    if (!is.null(fit_reg)) {
      is_model_valid <- TRUE
    } else {
      warning("fit_reg is NULL for loop ", loop_id)
    }
  } else {
    warning("Missing file for ", loop_id)
  }

  if (is_model_valid) {
    coef_mat <- fit_reg@coefficients
    rownames(coef_mat) <- colnames(fit_reg@data$X)
    test_mat <- fit_reg@test
    
    logLik_val <- fit_reg@logL
    BIC_val    <- fit_reg@BIC
    AIC_val    <- fit_reg@AIC
    iter_val   <- fit_reg@iter
    
    grad_mat <- fit_reg@gradient
    grad_norm_val <- sum(grad_mat^2, na.rm = TRUE)
    is_converged <- mean(grad_mat^2, na.rm = TRUE) <= 1e-04
    
    p_val_time <- get_coef_safe(test_mat, "time", "Pr(>wald)")
    if (is.na(p_val_time)) {
      fit_status <- "FAILED_PVAL_NA"
    } else if (!is_converged) {
      fit_status <- "NOT_CONVERGED"
    } else {
      fit_status <- "SUCCESS"
    }
    
  } else {
    coef_mat <- NULL
    test_mat <- NULL
    logLik_val = BIC_val = AIC_val = iter_val = NA_real_
    grad_norm_val = NA_real_
    is_converged = NA
    p_val_time = NA_real_
    fit_status <- "NO_MODEL_FILE"
  }

  row_dt <- data.table(
    loop_id = loop_id,
    
    # INTERCEPT
    intercept_alpha_x_out = get_coef_safe(coef_mat, "(Intercept)", "alpha_x_out"),
    intercept_alpha_x_A2  = get_coef_safe(coef_mat, "(Intercept)", "alpha_x_A2"),
    intercept_alpha_x_A1  = get_coef_safe(coef_mat, "(Intercept)", "alpha_x_A1"),
    intercept_beta_x_out  = get_coef_safe(coef_mat, "(Intercept)", "beta_x_out"),
    intercept_beta_x_A2   = get_coef_safe(coef_mat, "(Intercept)", "beta_x_A2"),
    intercept_beta_x_A1   = get_coef_safe(coef_mat, "(Intercept)", "beta_x_A1"),
  
    # TIME
    time_alpha_x_out = get_coef_safe(coef_mat, "time", "alpha_x_out"),
    time_alpha_x_A2  = get_coef_safe(coef_mat, "time", "alpha_x_A2"),
    time_alpha_x_A1  = get_coef_safe(coef_mat, "time", "alpha_x_A1"),
    time_beta_x_out  = get_coef_safe(coef_mat, "time", "beta_x_out"),
    time_beta_x_A2   = get_coef_safe(coef_mat, "time", "beta_x_A2"),
    time_beta_x_A1   = get_coef_safe(coef_mat, "time", "beta_x_A1"),
    
    intercept_wald   = get_coef_safe(test_mat, "(Intercept)", "wald value"),
    intercept_pval   = get_coef_safe(test_mat, "(Intercept)", "Pr(>wald)"),
    time_wald        = get_coef_safe(test_mat, "time", "wald value"),
    time_pval        = get_coef_safe(test_mat, "time", "Pr(>wald)"),
    
    converged     = is_converged,
    gradient_norm = grad_norm_val,
    fit_status     = fit_status,
  
    logLik      = logLik_val,
    BIC         = BIC_val,
    AIC         = AIC_val,
    iterations  = iter_val
  )

  results_list[[loop_id]] <- row_dt
}

results_dt <- rbindlist(results_list, fill = TRUE)
fwrite(results_dt, output_file, sep = "\t")
