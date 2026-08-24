library(data.table)

config_path <- file.path("config", "config.yml")
config <- yaml::yaml.load_file(config_path)

srcdir     <- config$paths$srcdir
resultsdir <- config$paths$resultsdir

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage: Rscript MGLMfit_GDM_check_and_corr.R <rho_value> <output_file_tsv_gz> <fit_dir> <synthetic_data_dir>")
}

rho_val        <- as.numeric(args[1])
output_file    <- args[2]
fit_dir        <- args[3]
synthetic_dir  <- args[4]

rho_str <- format(rho_val, nsmall = 1)
output_dir <- dirname(output_file)
if (!dir.exists(output_dir)) dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

N_CELLS <- c(500, 1000, 1500, 2000, 2500, 3000, 3500, 4000, 5000, 
             10000, 15000, 20000, 25000, 30000, 35000, 40000, 45000, 
             50000, 55000, 60000, 80000, 100000, 120000, 140000, 160000)
MUS <- c(1000, 2000, 3000, 4000, 5000, 6000)
SIZE_NEGBINOMS <- c("inf", "fixed", "0.5", "1.0", "1.5", "2.0", "2.5", "3.0", "3.5", "10")
ALPHA <- c(0.15)
BETA <- c(600)

n_sim <- 10000
set.seed(12345)

rGDM3 <- function(n, fit, sorted_names) {
  alpha <- sapply(sorted_names[1:2], function(nm) fit@estimate[paste0("alpha_", nm)])
  beta  <- sapply(sorted_names[1:2], function(nm) fit@estimate[paste0("beta_", nm)])

  p1 <- rbeta(n, alpha[1], beta[1])
  p2 <- rbeta(n, alpha[2], beta[2]) * (1 - p1)
  p3 <- pmax(0, 1 - p1 - p2)

  res <- cbind(p1, p2, p3)
  colnames(res) <- sorted_names
  res
}

safe_extract <- function(x, name) {
  if (is.null(x)) return(NA_real_)
  if (!(name %in% names(x))) return(NA_real_)
  as.numeric(x[name])
}

cor_results <- list()

total_n <- length(N_CELLS)
current_n_idx <- 0

cat("=== STARTING PROCESSING FOR RHO:", rho_str, "===\n")
flush.console()

for (n in N_CELLS) {
  current_n_idx <- current_n_idx + 1

  cat(sprintf("[%s] Rho: %s | Progress: %d/%d n_cells group (current n: %d)\n", 
              format(Sys.time(), "%H:%M:%S"), rho_str, current_n_idx, total_n, n))
  flush.console()
  
  for (mu in MUS) {
    cat(sprintf("  -> Processing mu: %d\n", mu))
    flush.console()
    
    for (size_nb in SIZE_NEGBINOMS) {
      for (alpha_val in ALPHA) {
        for (beta_val in BETA) {

          reads_filename <- paste0("synthetic_counts_n", n, "_mu", mu, "_sizeNB", size_nb, "_alpha", alpha_val, "_beta", beta_val, ".tsv.gz")
          reads_file     <- file.path(synthetic_dir, paste0("rho_", rho_str), reads_filename)
          
          fit_filename   <- paste0("fit_synthetic_n", n, "_mu", mu, "_sizeNB", size_nb, "_alpha", alpha_val, "_beta", beta_val, ".rds")
          fit_file       <- file.path(fit_dir, paste0("rho_", rho_str), fit_filename)

          base_dt <- data.table(
            rho = rho_val, n = n, mu = mu, size_nb = size_nb, alpha_param = alpha_val, beta_param = beta_val,
            spearman_rho = NA_real_, spearman_pvalue = NA_real_,
            pearson_r = NA_real_, pearson_pvalue = NA_real_,
            alpha_x_A1_est = NA_real_, alpha_x_A1_SE = NA_real_,
            alpha_x_A2_est = NA_real_, alpha_x_A2_SE = NA_real_,
            alpha_x_out_est = NA_real_, alpha_x_out_SE = NA_real_,
            beta_x_A1_est = NA_real_, beta_x_A1_SE = NA_real_,
            beta_x_A2_est = NA_real_, beta_x_A2_SE = NA_real_,
            beta_x_out_est = NA_real_, beta_x_out_SE = NA_real_,
            min_est_over_se = NA_real_
          )
          
          if (!file.exists(reads_file)) {
            cor_results[[length(cor_results) + 1]] <- base_dt
            next
          }
          
          fit <- tryCatch(readRDS(fit_file), error = function(e) NULL)
          if (is.null(fit) || is.null(fit@estimate)) {
            cor_results[[length(cor_results) + 1]] <- base_dt
            next
          }
          
          reads <- fread(reads_file, select = c("x_A1", "x_A2", "total_reads"))
          reads[, x_out := total_reads - x_A1 - x_A2]
          
          col_sums <- colSums(reads[, .(x_A1, x_A2, x_out)])
          ord <- order(col_sums, decreasing = TRUE)
          sorted_names <- names(col_sums)[ord]
          
          rm(reads)
          

          est <- fit@estimate
          se <- fit@SE
          
          alpha_est <- setNames(rep(NA_real_, 3), sorted_names)
          beta_est  <- setNames(rep(NA_real_, 3), sorted_names)
          alpha_se  <- setNames(rep(NA_real_, 3), sorted_names)
          beta_se   <- setNames(rep(NA_real_, 3), sorted_names)

          for (nm in sorted_names[1:2]) {
            alpha_est[nm] <- safe_extract(est, paste0("alpha_", nm))
            beta_est[nm]  <- safe_extract(est, paste0("beta_", nm))
          }
          
          alpha_se[sorted_names[1]] <- se[[1]]
          alpha_se[sorted_names[2]] <- se[[2]]
          beta_se[sorted_names[1]]  <- se[[3]]
          beta_se[sorted_names[2]]  <- se[[4]]
          
          est_vec <- c(alpha_est, beta_est)
          se_vec  <- c(alpha_se, beta_se)
          
          all_params_vec <- c(est_vec, se_vec)

          if (sum(!is.na(all_params_vec)) < 8) {
            base_dt[, `:=`(
              alpha_x_A1_est = alpha_est["x_A1"], alpha_x_A1_SE = alpha_se["x_A1"],
              alpha_x_A2_est = alpha_est["x_A2"], alpha_x_A2_SE = alpha_se["x_A2"],
              alpha_x_out_est = alpha_est["x_out"], alpha_x_out_SE = alpha_se["x_out"],
              beta_x_A1_est = beta_est["x_A1"], beta_x_A1_SE = beta_se["x_A1"],
              beta_x_A2_est = beta_est["x_A2"], beta_x_A2_SE = beta_se["x_A2"],
              beta_x_out_est = beta_est["x_out"], beta_x_out_SE = beta_se["x_out"]
            )]
            cor_results[[length(cor_results) + 1]] <- base_dt
            next
          }
          
          y_sim <- rGDM3(n_sim, fit, sorted_names)
          dt <- as.data.table(y_sim)
          
          ct_s <- suppressWarnings(cor.test(dt$x_A1, dt$x_A2, method = "spearman"))
          ct_p <- suppressWarnings(cor.test(dt$x_A1, dt$x_A2, method = "pearson"))
          
          ratio <- abs(est_vec) / se_vec
          ratio <- ratio[is.finite(ratio)]
          min_ratio_val <- ifelse(length(ratio) == 0, NA_real_, min(ratio))
          
          res_dt <- data.table(
            rho = rho_val, n = n, mu = mu, size_nb = size_nb, alpha_param = alpha_val, beta_param = beta_val,
            spearman_rho = as.numeric(ct_s$estimate),
            spearman_pvalue = ct_s$p.value,
            pearson_r = as.numeric(ct_p$estimate),
            pearson_pvalue = ct_p$p.value,
            alpha_x_A1_est = alpha_est["x_A1"], alpha_x_A1_SE = alpha_se["x_A1"],
            alpha_x_A2_est = alpha_est["x_A2"], alpha_x_A2_SE = alpha_se["x_A2"],
            alpha_x_out_est = alpha_est["x_out"], alpha_x_out_SE = alpha_se["x_out"],
            beta_x_A1_est = beta_est["x_A1"], beta_x_A1_SE = beta_se["x_A1"],
            beta_x_A2_est = beta_est["x_A2"], beta_x_A2_SE = beta_se["x_A2"],
            beta_x_out_est = beta_est["x_out"], beta_x_out_SE = beta_se["x_out"],
            min_est_over_se = min_ratio_val
          )
          
          cor_results[[length(cor_results) + 1]] <- res_dt
        }
      }
    }
  }
}

cat("=== WRITING RESULTS FOR RHO:", rho_str, "===\n")
flush.console()

final_dt <- rbindlist(cor_results)
fwrite(final_dt, output_file, sep = "\t", compress = "gzip")

cat("=== FINISHED PROCESSING FOR RHO:", rho_str, "===\n")
flush.console()
