library(data.table)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
  stop("Usage: Rscript cor_synthetic.R <rho_value> <init_value> <output_file_tsv_gz> <fit_base_dir> <synthetic_data_dir>")
}

rho_val        <- as.numeric(args[1])
init_val       <- args[2]
output_file    <- args[3]
fit_base_dir   <- args[4]
synthetic_dir  <- args[5]

rho_str <- format(rho_val, nsmall = 1)
source(file.path("src", "R", "MGLMfit_GDM_cor", "correlation_functions.R"))

fit_dir <- file.path(fit_base_dir, "MGLMfit_GDM", "synthetic_data", paste0("rho_", rho_str), paste0("init_", init_val))

N_CELLS <- c("500", "1000", "1500", "2000", "2500", "3000", "3500", "4000", "5000", 
             "10000", "15000", "20000", "25000", "30000", "35000", "40000", "45000", 
             "50000", "55000", "60000", "80000", "100000", "120000", "140000", "160000")
MUS <- c(1000, 2000, 3000, 4000, 5000, 6000)
SIZE_NEGBINOMS <- c("inf", "fixed", "0.5", "1.0", "1.5", "2.0", "2.5", "3.0", "3.5", "10")
ALPHA <- c(0.15)
BETA <- c(600)

cor_results <- list()

for (n in N_CELLS) {
  for (mu in MUS) {
    for (size_nb in SIZE_NEGBINOMS) {
      for (alpha_val in ALPHA) {
        for (beta_val in BETA) {
          
          reads_file <- file.path(synthetic_dir, paste0("rho_", rho_str), paste0("synthetic_counts_n", n, "_mu", mu, "_sizeNB", size_nb, "_alpha", alpha_val, "_beta", beta_val, ".tsv.gz"))
          fit_file   <- file.path(fit_dir, paste0("fit_synthetic_n", n, "_mu", mu, "_sizeNB", size_nb, "_alpha", alpha_val, "_beta", beta_val, ".rds"))
          
          meta <- list(rho = rho_val, n = n, mu = mu, size_nb = size_nb, alpha_param = alpha_val, beta_param = beta_val)
          
          if (!file.exists(reads_file)) {
            cor_results[[length(cor_results) + 1]] <- compute_gdm_correlations(NULL, NULL, metadata_list = meta)
            next
          }
          
          fit <- tryCatch(readRDS(fit_file), error = function(e) NULL)
          cor_results[[length(cor_results) + 1]] <- compute_gdm_correlations(fit, reads_file, mode = "all", metadata_list = meta)
        }
      }
    }
  }
}

fwrite(rbindlist(cor_results), output_file, sep = "\t", compress = "gzip")
