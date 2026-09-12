library(data.table)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
  stop("Usage: Rscript cor_synthetic_extended_plus.R <rho_value> <init_value> <seed_value> <output_file_tsv_gz> <fit_base_dir>")
}

rho_val      <- as.numeric(args[1])
init_val     <- args[2]
seed_val     <- as.integer(args[3])
output_file  <- args[4]
fit_base_dir <- args[5]

rho_str <- format(rho_val, nsmall = 1)
source(file.path("src", "R", "MGLMfit_GDM_cor", "correlation_functions.R"))

fit_dir <- file.path(fit_base_dir, "MGLMfit_GDM", "synthetic_data_extended_plus", paste0("rho_", rho_str), paste0("init_", init_val))

N_CELLS        <- c("1000", "2000", "5000", "10000", "20000")
MUS            <- c(1000, 2000, 3000, 4000, 5000)
SIZE_NEGBINOMS <- c("0.1", "0.2", "0.5", "1.0", "2.0", "inf", "fixed")

ALPHA <- 0.15
BETA  <- 600

cor_results <- list()

for (n in N_CELLS) {
  for (mu in MUS) {
    for (size_nb in SIZE_NEGBINOMS) {

      fit_file <- file.path(fit_dir, paste0("fit_synthetic_n", n, "_mu", mu, "_sizeNB", size_nb, "_alpha", ALPHA, "_beta", BETA, "_seed", seed_val, ".rds"))

      meta <- list(
        rho = rho_val, 
        seed_param = seed_val,
        n = as.integer(n), 
        mu = mu,
        size_nb = size_nb, 
        alpha_param = ALPHA, 
        beta_param = BETA
      )
      
      if (!file.exists(fit_file)) {
        cor_results[[length(cor_results) + 1]] <- compute_gdm_correlations_from_fit(NULL, metadata_list = meta)
        next
      }
      
      fit <- tryCatch(readRDS(fit_file), error = function(e) NULL)

      cor_results[[length(cor_results) + 1]] <- compute_gdm_correlations_from_fit(
        fit = fit, 
        metadata_list = meta
      )
    }
  }
}

out_dir <- dirname(output_file)
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
}

fwrite(rbindlist(cor_results), output_file, sep = "\t", compress = "gzip")
cat("Successfully computed extended correlations for rho:", rho_str, "seed:", seed_val, "\n")
