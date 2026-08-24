library(data.table)
library(MASS)

####################################################################################################
# Bivariate Beta-Binomial Generator
####################################################################################################
bivrbetabinom <- function(n, shape1, shape2, sizes, sigma = 0) {
    Sigma <- matrix(c(1, sigma, sigma, 1), ncol = 2)
    Z <- mvrnorm(n, rep(0, 2), Sigma)
    prob <- cbind(
        qbeta(pnorm(Z[, 1]), shape1 = shape1[1], shape2 = shape2[1]),
        qbeta(pnorm(Z[, 2]), shape1 = shape1[2], shape2 = shape2[2]))
    
    generated_norm_cor <- cor(Z[,1], Z[,2])
    prob_pearson_cor  <- cor(prob[,1], prob[,2], method = "pearson")
    prob_spearman_cor <- cor(prob[,1], prob[,2], method = "spearman")

    counts <- cbind(
        rbinom(n, size = sizes, prob = prob[, 1]),
        rbinom(n, size = sizes, prob = prob[, 2]))

    return(list(
        counts = counts,
        sigma = sigma,
        generated_norm_cor = generated_norm_cor,
        prob_pearson_cor = prob_pearson_cor,
        prob_spearman_cor = prob_spearman_cor
    ))
}

####################################################################################################

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Missing arguments!\nUsage: Rscript generate_synthetic_data_extended.R <rho> <seed>")
}

rho_val   <- as.numeric(args[1])
seed_val  <- as.integer(args[2])

#n_cells_list <- c(
#  500, 1000, 1500, 2000, 2500, 3000, 3500, 4000, 5000,
#  10000, 20000, 40000, 60000, 80000, 100000, 120000, 140000, 160000
#)
#mus_list        <- c(1000, 2000, 3000, 4000, 5000, 6000)
#size_nb_strings <- c("inf", "fixed", "0.1", "0.2", "0.3", "0.4", "0.5", "0.6", "0.7", "0.8", "0.9", "1.0", "1.5", "2.0", "2.5", "3.0", "3.5", "10.0")

n_cells_list <- c(1000, 2000, 5000, 10000, 20000)
mus_list        <- c(5000)
size_nb_strings <- c("inf", "fixed", "0.5", "1.0", "2.0")

alpha_val <- 0.15
beta_val  <- 600

config_path <- file.path("config", "config.yml")
config      <- yaml::yaml.load_file(config_path)
resultsdir  <- config$paths$resultsdir

rho_str    <- format(rho_val, nsmall = 1)
output_dir <- file.path(resultsdir, "synthetic_data_extended", paste0("rho_", rho_str))
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
}

set.seed(42 + seed_val)

for (n in n_cells_list) {
  for (mu_val in mus_list) {
    for (size_nb_str in size_nb_strings) {

      if (size_nb_str == "fixed") {
        sizes <- rep(mu_val, n)
      } else if (size_nb_str == "inf") {
        sizes <- rpois(n, lambda = mu_val)
      } else {
        size_nb_val <- as.numeric(size_nb_str)
        if (is.na(size_nb_val)) next
        sizes <- rnbinom(n, size = size_nb_val, mu = mu_val)
      }
      
      sizes <- pmax(sizes, 1)
      
      sim <- bivrbetabinom(
        n = n,
        shape1 = c(alpha_val, alpha_val),
        shape2 = c(beta_val, beta_val),
        sizes = sizes,
        sigma = rho_val
      )
      
      counts   <- sim$counts
      x_A1_syn <- counts[, 1]
      x_A2_syn <- counts[, 2]
      
      if (any(x_A1_syn > sizes) || any(x_A2_syn > sizes)) {
        stop("Generated counts violate constraint: x_A1 or x_A2 > total_reads")
      }
      
      out <- data.table(
        total_reads = sizes,
        x_A1 = x_A1_syn,
        x_A2 = x_A2_syn
      )
      
      out[, sigma := NA_real_]
      out[, generated_norm_cor := NA_real_]
      out[, prob_pearson_cor := NA_real_]
      out[, prob_spearman_cor := NA_real_]
      
      out[1, `:=`(
        sigma = sim$sigma,
        generated_norm_cor = sim$generated_norm_cor,
        prob_pearson_cor = sim$prob_pearson_cor,
        prob_spearman_cor = sim$prob_spearman_cor
      )]
      
      filename <- paste0("synthetic_counts_n", n, "_mu", mu_val, "_sizeNB", size_nb_str, "_alpha", alpha_val, "_beta", beta_val, "_seed", seed_val, ".tsv.gz")
      output_file <- file.path(output_dir, filename)
      
      fwrite(out, output_file, sep = "\t")
    }
  }
}

done_file <- file.path(output_dir, paste0("sim_seed", seed_val, ".done"))
file.create(done_file)
cat("Successfully finished batch for rho:", rho_str, "and seed:", seed_val, "\n")
