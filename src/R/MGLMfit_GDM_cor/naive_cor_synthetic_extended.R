library(data.table)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage: Rscript naive_cor_synthetic_extended.R <rho_value> <seed_value> <output_file_tsv_gz> <synthetic_data_dir>")
}

rho_val       <- as.numeric(args[1])
seed_val      <- as.integer(args[2])
output_file   <- args[3]
synthetic_dir <- args[4]

rho_str <- format(rho_val, nsmall = 1)

N_CELLS <- c(
  "500", "1000", "1500", "2000", "2500", "3000", "3500", "4000", "5000",
  "10000", "20000", "40000", "60000", "80000", "100000", "120000", "140000", "160000"
)
MUS <- c(1000, 2000, 3000, 4000, 5000, 6000)
SIZE_NEGBINOMS <- c("inf", "fixed", "0.1", "0.2", "0.3", "0.4", "0.5", "0.6", "0.7", "0.8", "0.9", "1.0", "1.5", "2.0", "2.5", "3.0", "3.5", "10.0")

#N_CELLS <- c("1000", "2000", "5000", "10000", "20000")
#MUS <- c(5000)
#SIZE_NEGBINOMS <- c("inf", "fixed", "0.5", "1.0", "2.0")

ALPHA <- c(0.15)
BETA <- c(600)

cor_results <- list()

for (n in N_CELLS) {
  for (mu in MUS) {
    for (size_nb in SIZE_NEGBINOMS) {
      for (alpha_val in ALPHA) {
        for (beta_val in BETA) {
          
          filename <- paste0("synthetic_counts_n", n, "_mu", mu, "_sizeNB", size_nb, "_alpha", alpha_val, "_beta", beta_val, "_seed", seed_val, ".tsv.gz")
          reads_file <- file.path(synthetic_dir, paste0("rho_", rho_str), filename)
          
          meta <- list(rho = rho_val, seed = seed_val, n = n, mu = mu, size_nb = size_nb, alpha_param = alpha_val, beta_param = beta_val)
          base_dt <- as.data.table(meta)
          base_dt[, `:=`(
            spearman_rho = NA_real_, spearman_pvalue = NA_real_,
            pearson_r = NA_real_, pearson_pvalue = NA_real_
          )]
          
          if (file.exists(reads_file)) {
            reads <- fread(reads_file, select = c("x_A1", "x_A2"))
            if (nrow(reads) >= 3) {
              ct_s <- suppressWarnings(tryCatch(cor.test(reads$x_A1, reads$x_A2, method = "spearman"), error = function(e) NULL))
              ct_p <- suppressWarnings(tryCatch(cor.test(reads$x_A1, reads$x_A2, method = "pearson"), error = function(e) NULL))
              
              if (!is.null(ct_s) && !is.null(ct_p)) {
                base_dt[, `:=`(
                  spearman_rho = as.numeric(ct_s$estimate),
                  spearman_pvalue = ct_s$p.value,
                  pearson_r = as.numeric(ct_p$estimate),
                  pearson_pvalue = ct_p$p.value
                )]
              }
            }
          }
          cor_results[[length(cor_results) + 1]] <- base_dt
        }
      }
    }
  }
}

output_dir <- dirname(output_file)
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
}

fwrite(rbindlist(cor_results), output_file, sep = "\t", compress = "gzip")
cat("Finished naive correlation calculation for rho:", rho_str, "seed:", seed_val, "\n")
