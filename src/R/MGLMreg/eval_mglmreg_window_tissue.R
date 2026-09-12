library(data.table)
library(MGLM)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript eval_mglmreg_window_tissue.R <tw> <output_file> [init_type]")
}

tw          <- args[1]
output_file <- args[2]
init_type   <- if (length(args) >= 3) args[3] else "default"
fit_init_type <- "1e-6"

message("==================================================")
message(sprintf("Start przetwarzania okna czasowego: %s (10h+ tissue)", tw))
message(sprintf("Typ inicjalizacji: %s", init_type))
message(sprintf("Plik docelowy: %s", output_file))
message("==================================================")

time_set   <- "10h+"
model_type <- "time_tissue"

config_path <- file.path("config", "config.yml")
config      <- yaml::yaml.load_file(config_path)

resultsdir <- config$paths$resultsdir
srcdir     <- config$paths$srcdir

# Ładowanie lokalnego pakietu MGLM
pkg_path <- file.path(srcdir, "MGLM")
if (dir.exists(pkg_path)) {
  r_files <- list.files(file.path(pkg_path, "R"), pattern = "\\.R$", full.names = TRUE)
  for (f in r_files) source(f)
}

tw_bounds <- as.numeric(strsplit(tw, "-")[[1]])
mid_time  <- mean(tw_bounds)
message(sprintf("Środkowy punkt czasowy (mid_time): %.2f h", mid_time))

get_coef_safe <- function(mat, row_name, col_name) {
  if (is.null(mat) || (!is.matrix(mat) && !is.data.frame(mat))) {
    return(NA_real_)
  }
  if (row_name %in% rownames(mat) && col_name %in% colnames(mat)) {
    val <- mat[row_name, col_name]
    if (length(val) == 1 && !is.na(val)) return(as.numeric(val))
  }
  r_idx <- which(rownames(mat) == row_name)
  if (length(r_idx) > 0) {
    c_idx <- grep(col_name, colnames(mat), ignore.case = TRUE)
    if (length(c_idx) > 0) {
      val <- mat[r_idx[1], c_idx[1]]
      if (length(val) == 1 && !is.na(val)) return(as.numeric(val))
    }
  }
  return(NA_real_)
}

sim_gdm_rho <- function(sorted_names, alpha_vec, beta_vec, n_sim = 10000) {
  if (any(is.na(alpha_vec)) || any(is.na(beta_vec)) || any(alpha_vec <= 0) || any(beta_vec <= 0)) {
    return(list(rho = NA_real_, pval = NA_real_))
  }
  
  p1 <- rbeta(n_sim, alpha_vec[1], beta_vec[1])
  p2 <- rbeta(n_sim, alpha_vec[2], beta_vec[2]) * (1 - p1)
  p3 <- pmax(0, 1 - p1 - p2)
  
  dt_sim <- data.table(p1, p2, p3)
  setnames(dt_sim, sorted_names)
  
  if (!all(c("x_A1", "x_A2") %in% names(dt_sim))) {
    return(list(rho = NA_real_, pval = NA_real_))
  }
  
  ct <- suppressWarnings(cor.test(dt_sim$x_A1, dt_sim$x_A2, method = "spearman"))
  return(list(rho = as.numeric(ct$estimate), pval = ct$p.value))
}

# Plik korelacji dla poszczególnych populacji
fit_cor_file <- file.path(resultsdir, "MGLMfit_GDM_cor", "real_data", "pops", paste0("init_", fit_init_type), paste0("cor_", tw, ".tsv.gz"))
if (!file.exists(fit_cor_file)) {
  stop("Nie znaleziono pliku MGLMfit cor populacji: ", fit_cor_file)
}

dt_fit <- fread(fit_cor_file)

req_cols <- c("time_window", "population", "loop_id", "min_est_over_se", 
              "spearman_rho", "spearman_pvalue",
              "alpha_x_A1_est", "alpha_x_A2_est", "alpha_x_out_est",
              "beta_x_A1_est", "beta_x_A2_est", "beta_x_out_est")

existing_cols <- intersect(req_cols, names(dt_fit))
dt_res <- dt_fit[, ..existing_cols]

setnames(dt_res, 
         old = c("min_est_over_se", "spearman_rho", "spearman_pvalue",
                 "alpha_x_A1_est", "alpha_x_A2_est", "alpha_x_out_est",
                 "beta_x_A1_est", "beta_x_A2_est", "beta_x_out_est"),
         new = c("fit_min_est_over_se", "fit_spearman_rho", "fit_spearman_pvalue",
                 "fit_alpha_x_A1_est", "fit_alpha_x_A2_est", "fit_alpha_x_out_est",
                 "fit_beta_x_A1_est", "fit_beta_x_A2_est", "fit_beta_x_out_est"),
         skip_absent = TRUE)

reg_dir <- file.path(resultsdir, "MGLMreg", sprintf("GDM_%s_%s_%s", time_set, model_type, init_type))
message(sprintf("Katalog z plikami RDS MGLMreg: %s", reg_dir))

dt_res[, `:=`(
  tissue              = NA_character_,
  reg_status          = NA_character_,
  reg_alpha_x_A1_est  = NA_real_, reg_beta_x_A1_est  = NA_real_,
  reg_alpha_x_A2_est  = NA_real_, reg_beta_x_A2_est  = NA_real_,
  reg_alpha_x_out_est = NA_real_, reg_beta_x_out_est = NA_real_,
  reg_pred_p_x_A1     = NA_real_, reg_pred_p_x_A2     = NA_real_, reg_pred_p_x_out = NA_real_,
  reg_spearman_rho    = NA_real_, reg_spearman_pvalue = NA_real_
)]

# Ustalenie nazwy tkanki dokładnie tak jak w GDM_unified_reg.R (sub("^[^_]+_", "", population))
dt_res[, tissue := sub("^[^_]+_", "", population)]

unique_loops <- unique(dt_res$loop_id)
n_loops      <- length(unique_loops)
message(sprintf("Liczba unikalnych pętli do przetworzenia: %d (łączna liczba wierszy: %d)", n_loops, nrow(dt_res)))

for (l_idx in seq_along(unique_loops)) {
  loop    <- unique_loops[l_idx]
  verbose <- (l_idx <= 3)
  
  if (verbose) {
    message(sprintf("\n--- [DEBUG LOOP %d/%d] ID: %s ---", l_idx, n_loops, loop))
  } else if (l_idx %% 50 == 0) {
    message(sprintf("[%s] Okno: %s | Pętla [%d/%d]: %s", format(Sys.time(), "%H:%M:%S"), tw, l_idx, n_loops, loop))
  }
  flush.console()
  
  row_indices <- which(dt_res$loop_id == loop)
  reg_file    <- file.path(reg_dir, paste0("fit_reg_", loop, ".rds"))
  
  if (!file.exists(reg_file)) {
    dt_res[row_indices, reg_status := "NO_MODEL_FILE"]
    if (verbose) message("  [!] Status: NO_MODEL_FILE (brak pliku)")
    next
  }
  
  fit_reg <- tryCatch(readRDS(reg_file), error = function(e) NULL)
  
  if (is.null(fit_reg) || !inherits(fit_reg, "MGLMreg")) {
    dt_res[row_indices, reg_status := "NO_MODEL_FILE"]
    if (verbose) message("  [!] Status: NO_MODEL_FILE (błąd odczytu RDS / obiekt NULL)")
    next
  }
  
  test_mat <- fit_reg@test
  grad_mat <- fit_reg@gradient
  
  if (is.null(grad_mat)) {
    is_converged <- FALSE
  } else {
    is_converged <- mean(grad_mat^2, na.rm = TRUE) <= 1e-04
  }
  
  p_val_time <- get_coef_safe(test_mat, "time", "Pr(>wald)")
  
  if (is.na(p_val_time)) {
    fit_status <- "FAILED_PVAL_NA"
  } else if (!is_converged) {
    fit_status <- "NOT_CONVERGED"
  } else {
    fit_status <- "SUCCESS"
  }
  
  dt_res[row_indices, reg_status := fit_status]
  
  if (verbose) {
    message(sprintf("  [+] Status: %s (is_converged=%s, p_val_time=%s)", 
                    fit_status, is_converged, p_val_time))
  }
  
  if (fit_status != "SUCCESS") next
  
  coef_mat <- fit_reg@coefficients
  if (is.null(coef_mat)) next
  
  if (!is.null(test_mat) && !is.null(rownames(test_mat))) {
    rownames(coef_mat) <- rownames(test_mat)
  }
  
  coef_cols <- colnames(coef_mat)
  coef_rows <- rownames(coef_mat)
  
  for (idx in row_indices) {
    t_name <- dt_res$tissue[idx]
    tissue_row_name <- paste0("tissue", t_name)
    
    reg_alphas <- list()
    reg_betas  <- list()
    
    for (cat_name in c("x_out", "x_A2", "x_A1")) {
      col_a <- paste0("alpha_", cat_name)
      col_b <- paste0("beta_", cat_name)
      
      if (col_a %in% coef_cols) {
        eta_a <- 0
        if ("(Intercept)" %in% coef_rows) eta_a <- eta_a + coef_mat["(Intercept)", col_a]
        if ("time" %in% coef_rows)        eta_a <- eta_a + mid_time * coef_mat["time", col_a]
        if (tissue_row_name %in% coef_rows) eta_a <- eta_a + coef_mat[tissue_row_name, col_a]
        reg_alphas[[cat_name]] <- exp(eta_a)
      } else {
        reg_alphas[[cat_name]] <- NA_real_
      }
      
      if (col_b %in% coef_cols) {
        eta_b <- 0
        if ("(Intercept)" %in% coef_rows) eta_b <- eta_b + coef_mat["(Intercept)", col_b]
        if ("time" %in% coef_rows)        eta_b <- eta_b + mid_time * coef_mat["time", col_b]
        if (tissue_row_name %in% coef_rows) eta_b <- eta_b + coef_mat[tissue_row_name, col_b]
        reg_betas[[cat_name]] <- exp(eta_b)
      } else {
        reg_betas[[cat_name]] <- NA_real_
      }
    }
    
    dt_res[idx, `:=`(
      reg_alpha_x_A1_est  = reg_alphas[["x_A1"]], reg_beta_x_A1_est  = reg_betas[["x_A1"]],
      reg_alpha_x_A2_est  = reg_alphas[["x_A2"]], reg_beta_x_A2_est  = reg_betas[["x_A2"]],
      reg_alpha_x_out_est = reg_alphas[["x_out"]], reg_beta_x_out_est = reg_betas[["x_out"]]
    )]
    
    # Budowanie macierzy dla funkcji predict()
    x_vec <- numeric(length(coef_rows))
    names(x_vec) <- coef_rows
    if ("(Intercept)" %in% names(x_vec))   x_vec["(Intercept)"] <- 1
    if ("time" %in% names(x_vec))          x_vec["time"] <- mid_time
    if (tissue_row_name %in% names(x_vec)) x_vec[tissue_row_name] <- 1
    
    newdata_mat <- matrix(x_vec, nrow = 1)
    
    pred_prob <- tryCatch(predict(fit_reg, newdata = newdata_mat), error = function(e) {
      if (verbose) message(sprintf("  [!] BŁĄD W predict(): %s", e$message))
      return(NULL)
    })
    
    if (!is.null(pred_prob) && is.matrix(pred_prob)) {
      if ("x_A1" %in% colnames(pred_prob))  dt_res[idx, reg_pred_p_x_A1  := pred_prob[1, "x_A1"]]
      if ("x_A2" %in% colnames(pred_prob))  dt_res[idx, reg_pred_p_x_A2  := pred_prob[1, "x_A2"]]
      if ("x_out" %in% colnames(pred_prob)) dt_res[idx, reg_pred_p_x_out := pred_prob[1, "x_out"]]
    }
    
    alpha_cols <- grep("^alpha_", coef_cols, value = TRUE)
    step_cats  <- sub("^alpha_", "", alpha_cols)
    
    if ("x_out" %in% step_cats) {
      step2_cat    <- setdiff(step_cats, "x_out")
      leftover_cat <- setdiff(c("x_out", "x_A1", "x_A2"), step_cats)
      
      sorted_names <- c("x_out", step2_cat, leftover_cat)
      
      a_vec <- c(reg_alphas[["x_out"]], reg_alphas[[step2_cat]])
      b_vec <- c(reg_betas[["x_out"]],  reg_betas[[step2_cat]])
      
      sim_res <- sim_gdm_rho(sorted_names, a_vec, b_vec, n_sim = 10000)
      
      dt_res[idx, `:=`(
        reg_spearman_rho    = sim_res$rho,
        reg_spearman_pvalue = sim_res$pval
      )]
    }
  }
}

if ("time_window" %in% names(dt_res)) {
  setcolorder(dt_res, c("time_window", "population", "tissue", "loop_id", "reg_status"))
}

out_dir <- dirname(output_file)
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
}

fwrite(dt_res, output_file, sep = "\t", compress = "gzip")
message(sprintf("Zakończono! Zapisano wyniki do: %s", output_file))
