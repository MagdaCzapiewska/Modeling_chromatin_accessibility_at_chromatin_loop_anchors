library(data.table)

args <- commandArgs(trailingOnly = TRUE)
#if (length(args) < 6) {
if (length(args) < 5) {
  #stop("Usage: Rscript cor_synthetic_extended.R <rho_value> <init_value> <seed_value> <output_file_tsv_gz> <fit_base_dir> <synthetic_data_dir>")
  stop("Usage: Rscript cor_synthetic_extended.R <rho_value> <init_value> <seed_value> <output_file_tsv_gz> <fit_base_dir>")
}

rho_val        <- as.numeric(args[1])
init_val       <- args[2]
seed_val       <- as.integer(args[3])
output_file    <- args[4]
fit_base_dir   <- args[5]
#synthetic_dir  <- args[6]

rho_str <- format(rho_val, nsmall = 1)
source(file.path("src", "R", "MGLMfit_GDM_cor", "correlation_functions.R"))

fit_dir <- file.path(fit_base_dir, "MGLMfit_GDM", "synthetic_data_extended", paste0("rho_", rho_str), paste0("init_", init_val))

#N_CELLS <- c(
#  "500", "1000", "1500", "2000", "2500", "3000", "3500", "4000", "5000",
#  "10000", "20000", "40000", "60000", "80000", "100000", "120000", "140000", "160000"
#)
#MUS <- c(1000, 2000, 3000, 4000, 5000, 6000)
#SIZE_NEGBINOMS <- c("inf", "fixed", "0.1", "0.2", "0.3", "0.4", "0.5", "0.6", "0.7", "0.8", "0.9", "1.0", "1.5", "2.0", "2.5", "3.0", "3.5", "10.0")

N_CELLS <- c("1000", "2000", "5000", "10000", "20000")
MUS <- c(5000)
SIZE_NEGBINOMS <- c("inf", "fixed", "0.5", "1.0", "2.0")

ALPHA <- 0.15
BETA <- 600

cor_results <- list()

for (n in N_CELLS) {
  for (mu in MUS) {
    for (size_nb in SIZE_NEGBINOMS) {

      #file_basename <- paste0("synthetic_counts_n", n, "_mu", mu, "_sizeNB", size_nb, "_alpha", ALPHA, "_beta", BETA, "_seed", seed_val, ".tsv.gz")
      #reads_file <- file.path(synthetic_dir, paste0("rho_", rho_str), file_basename)
      
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
      
      #if (!file.exists(reads_file)) {
      if (!file.exists(fit_file)) {
        #cor_results[[length(cor_results) + 1]] <- compute_gdm_correlations(NULL, NULL, metadata_list = meta)
        cor_results[[length(cor_results) + 1]] <- compute_gdm_correlations_from_fit(NULL, metadata_list = meta)
        next
      }
      
      fit <- tryCatch(readRDS(fit_file), error = function(e) NULL)

      #cor_results[[length(cor_results) + 1]] <- compute_gdm_correlations(
      #  fit = fit, 
      #  reads_file = reads_file, 
      #  mode = "all", 
      #  metadata_list = meta
      #)
      cor_results[[length(cor_results) + 1]] <- compute_gdm_correlations_from_fit(
        fit = fit, 
        metadata_list = meta
      )
    }
  }
}

fwrite(rbindlist(cor_results), output_file, sep = "\t", compress = "gzip")
cat("Successfully computed extended correlations for rho:", rho_str, "seed:", seed_val, "\n")
