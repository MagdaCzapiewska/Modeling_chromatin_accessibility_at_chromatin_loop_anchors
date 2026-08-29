library(data.table)
library(MGLM)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript eval_mglmreg_window.R <tw> <output_file>")
}

tw          <- args[1]
output_file <- args[2]

message("==================================================")
message(sprintf("Start przetwarzania okna czasowego: %s", tw))
message(sprintf("Plik docelowy: %s", output_file))
message("==================================================")

# Sztywno ustawione warianty
time_set   <- "0h+"
model_type <- "time"
init_type  <- "default"

config_path <- file.path("config", "config.yml")
config      <- yaml::yaml.load_file(config_path)

resultsdir <- config$paths$resultsdir
srcdir     <- config$paths$srcdir

# Ładowanie pakietu MGLM ze źródeł lokalnych (jeśli istnieją)
pkg_path <- file.path(srcdir, "MGLM")
if (dir.exists(pkg_path)) {
  r_files <- list.files(file.path(pkg_path, "R"), pattern = "\\.R$", full.names = TRUE)
  for (f in r_files) source(f)
}

# Obliczenie środkowego punktu czasowego (mid_time) dla danego okna
tw_bounds <- as.numeric(strsplit(tw, "-")[[1]])
mid_time  <- mean(tw_bounds)
message(sprintf("Środkowy punkt czasowy (mid_time): %.2f h", mid_time))

# Bezpieczne wyciąganie wartości z macierzy/tabeli testowej
get_coef_safe <- function(mat, row_name, col_name) {
  if (is.null(mat) || (!is.matrix(mat) && !is.data.frame(mat))) {
    return(NA_real_)
  }
  if (row_name %in% rownames(mat) && col_name %in% colnames(mat)) {
    val <- mat[row_name, col_name]
    if (length(val) == 1 && !is.na(val)) return(as.numeric(val))
  }
  # Dopasowanie elastyczne (np. brak wielkości liter)
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

# Funkcja do generowania korelacji ze zmapowanych parametrów GDM
sim_gdm_rho <- function(sorted_names, alpha_vec, beta_vec, n_sim = 10000) {
  if (any(is.na(alpha_vec)) || any(is.na(beta_vec)) || any(alpha_vec <= 0) || any(beta_vec <= 0)) {
    return(list(rho = NA_real_, pval = NA_real_))
  }
  
  # Krok 1: zawsze x_out
  p1 <- rbeta(n_sim, alpha_vec[1], beta_vec[1])
  # Krok 2: druga kategoria w kaskadzie (x_A1 lub x_A2)
  p2 <- rbeta(n_sim, alpha_vec[2], beta_vec[2]) * (1 - p1)
  # Krok 3: resztkowa kategoria bazowa
  p3 <- pmax(0, 1 - p1 - p2)
  
  dt_sim <- data.table(p1, p2, p3)
  setnames(dt_sim, sorted_names)
  
  if (!all(c("x_A1", "x_A2") %in% names(dt_sim))) {
    return(list(rho = NA_real_, pval = NA_real_))
  }
  
  ct <- suppressWarnings(cor.test(dt_sim$x_A1, dt_sim$x_A2, method = "spearman"))
  return(list(rho = as.numeric(ct$estimate), pval = ct$p.value))
}

# Wczytanie wyników MGLMfit
fit_cor_file <- file.path(resultsdir, "MGLMfit_GDM_cor", "real_data", "all", "init_1e-6", paste0("cor_", tw, ".tsv.gz"))
if (!file.exists(fit_cor_file)) {
  stop("Nie znaleziono pliku MGLMfit cor: ", fit_cor_file)
}

dt_fit <- fread(fit_cor_file)

req_cols <- c("time_window", "loop_id", "min_est_over_se", 
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

reg_dir <- file.path(resultsdir, "MGLMreg", "GDM_0h+_time_default")

dt_res[, `:=`(
  reg_status          = NA_character_,
  reg_alpha_x_A1_est  = NA_real_, reg_beta_x_A1_est  = NA_real_,
  reg_alpha_x_A2_est  = NA_real_, reg_beta_x_A2_est  = NA_real_,
  reg_alpha_x_out_est = NA_real_, reg_beta_x_out_est = NA_real_,
  reg_pred_p_x_A1     = NA_real_, reg_pred_p_x_A2     = NA_real_, reg_pred_p_x_out = NA_real_,
  reg_spearman_rho    = NA_real_, reg_spearman_pvalue = NA_real_
)]

n_rows <- nrow(dt_res)

for (i in seq_len(n_rows)) {
  loop    <- dt_res$loop_id[i]
  verbose <- (i <= 3)
  
  if (verbose) {
    message(sprintf("\n--- [DEBUG LOOP %d/%d] ID: %s ---", i, n_rows, loop))
  } else if (i %% 100 == 0) {
    message(sprintf("[%s] Okno: %s | Pętla [%d/%d]: %s", format(Sys.time(), "%H:%M:%S"), tw, i, n_rows, loop))
  }
  flush.console()
  
  reg_file <- file.path(reg_dir, paste0("fit_reg_", loop, ".rds"))
  
  # Brak pliku modelowego
  if (!file.exists(reg_file)) {
    dt_res[i, reg_status := "NO_MODEL_FILE"]
    if (verbose) message("  [!] Status: NO_MODEL_FILE (brak pliku)")
    next
  }
  
  fit_reg <- tryCatch(readRDS(reg_file), error = function(e) NULL)
  
  # Brak poprawnego obiektu po odczycie
  if (is.null(fit_reg) || !inherits(fit_reg, "MGLMreg")) {
    dt_res[i, reg_status := "NO_MODEL_FILE"]
    if (verbose) message("  [!] Status: NO_MODEL_FILE (błąd odczytu RDS / obiekt NULL)")
    next
  }
  
  test_mat <- fit_reg@test
  grad_mat <- fit_reg@gradient
  
  # Obliczenie wskaźników zbieżności i testów Walda
  if (is.null(grad_mat)) {
    is_converged <- FALSE
  } else {
    grad_norm_val <- sum(grad_mat^2, na.rm = TRUE)
    is_converged  <- mean(grad_mat^2, na.rm = TRUE) <= 1e-04
  }
  
  p_val_time <- get_coef_safe(test_mat, "time", "Pr(>wald)")
  
  if (is.na(p_val_time)) {
    fit_status <- "FAILED_PVAL_NA"
  } else if (!is_converged) {
    fit_status <- "NOT_CONVERGED"
  } else {
    fit_status <- "SUCCESS"
  }
  
  dt_res[i, reg_status := fit_status]
  
  if (verbose) {
    message(sprintf("  [+] Status: %s (is_converged=%s, p_val_time=%s)", 
                    fit_status, is_converged, p_val_time))
  }
  
  # Jeśli model nie odniósł sukcesu, pomijamy wyliczanie wartości
  #if (fit_status != "SUCCESS") next
  
  coef_mat <- fit_reg@coefficients
  if (is.null(coef_mat)) next
  
  # Ustawienie nazw wierszy z test_mat
  if (!is.null(test_mat) && !is.null(rownames(test_mat))) {
    rownames(coef_mat) <- rownames(test_mat)
  }
  
  coef_cols <- colnames(coef_mat)
  coef_rows <- rownames(coef_mat)
  
  reg_alphas <- list()
  reg_betas  <- list()
  
  # Wyznaczenie estymat parametrów w czasie mid_time: exp(intercept + time * mid_time)
  for (cat_name in c("x_out", "x_A2", "x_A1")) {
    col_a <- paste0("alpha_", cat_name)
    col_b <- paste0("beta_", cat_name)
    
    if (col_a %in% coef_cols && "(Intercept)" %in% coef_rows && "time" %in% coef_rows) {
      eta_a <- coef_mat["(Intercept)", col_a] + mid_time * coef_mat["time", col_a]
      reg_alphas[[cat_name]] <- exp(eta_a)
    } else {
      reg_alphas[[cat_name]] <- NA_real_
    }
    
    if (col_b %in% coef_cols && "(Intercept)" %in% coef_rows && "time" %in% coef_rows) {
      eta_b <- coef_mat["(Intercept)", col_b] + mid_time * coef_mat["time", col_b]
      reg_betas[[cat_name]] <- exp(eta_b)
    } else {
      reg_betas[[cat_name]] <- NA_real_
    }
  }
  
  dt_res[i, `:=`(
    reg_alpha_x_A1_est  = reg_alphas[["x_A1"]], reg_beta_x_A1_est  = reg_betas[["x_A1"]],
    reg_alpha_x_A2_est  = reg_alphas[["x_A2"]], reg_beta_x_A2_est  = reg_betas[["x_A2"]],
    reg_alpha_x_out_est = reg_alphas[["x_out"]], reg_beta_x_out_est = reg_betas[["x_out"]]
  )]
  
  # Predict – macierz numeryczna z odpowiednimi wymiarami
  if (!is.null(coef_rows) && "(Intercept)" %in% coef_rows && "time" %in% coef_rows) {
    newdata_mat <- matrix(c(1, mid_time), nrow = 1)
  } else {
    newdata_mat <- matrix(mid_time, nrow = 1)
  }
  
  pred_prob <- tryCatch(predict(fit_reg, newdata = newdata_mat), error = function(e) {
    if (verbose) message(sprintf("  [!] BŁĄD W predict(): %s", e$message))
    return(NULL)
  })
  
  if (!is.null(pred_prob) && is.matrix(pred_prob)) {
    if ("x_A1" %in% colnames(pred_prob))  dt_res[i, reg_pred_p_x_A1  := pred_prob[1, "x_A1"]]
    if ("x_A2" %in% colnames(pred_prob))  dt_res[i, reg_pred_p_x_A2  := pred_prob[1, "x_A2"]]
    if ("x_out" %in% colnames(pred_prob)) dt_res[i, reg_pred_p_x_out := pred_prob[1, "x_out"]]
  }
  
  # Wyznaczenie porządku kaskady GDM z gwarancją x_out na 1. pozycji
  alpha_cols <- grep("^alpha_", coef_cols, value = TRUE)
  step_cats  <- sub("^alpha_", "", alpha_cols)
  
  if ("x_out" %in% step_cats) {
    step2_cat    <- setdiff(step_cats, "x_out")
    leftover_cat <- setdiff(c("x_out", "x_A1", "x_A2"), step_cats)
    
    sorted_names <- c("x_out", step2_cat, leftover_cat)
    
    a_vec <- c(reg_alphas[["x_out"]], reg_alphas[[step2_cat]])
    b_vec <- c(reg_betas[["x_out"]],  reg_betas[[step2_cat]])
    
    sim_res <- sim_gdm_rho(sorted_names, a_vec, b_vec, n_sim = 10000)
    
    dt_res[i, `:=`(
      reg_spearman_rho    = sim_res$rho,
      reg_spearman_pvalue = sim_res$pval
    )]
  }
}

if ("time_window" %in% names(dt_res)) {
  setcolorder(dt_res, c("time_window", "loop_id", "reg_status"))
}

out_dir <- dirname(output_file)
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
}

fwrite(dt_res, output_file, sep = "\t", compress = "gzip")
message(sprintf("Zakończono! Zapisano wyniki do: %s", output_file))
