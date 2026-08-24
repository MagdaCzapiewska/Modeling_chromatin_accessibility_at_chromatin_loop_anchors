library(data.table)
library(MASS)

####################################################################################################

# Random generation of ‘n’ observations from two beta-binomial mixture
# distributions with parameters ‘shape1’, ‘shape2’ and ‘sizes’, with the
# covariation ‘sigma’ used in the construction of bivariate distribution.

bivrbetabinom <- function(n, shape1, shape2, sizes, sigma = 0) {
    Sigma <- matrix(c(1, sigma, sigma, 1), ncol = 2)
    Z <- mvrnorm(n, rep(0, 2), Sigma)
    prob <- cbind(
        qbeta(pnorm(Z[, 1]), shape1 = shape1[1], shape2 = shape2[1]),
        qbeta(pnorm(Z[, 2]), shape1 = shape1[2], shape2 = shape2[2]))
    
    generated_norm_cor <- cor(Z[,1], Z[,2])
    prob_pearson_cor  <- cor(prob[,1], prob[,2], method = "pearson")
    prob_spearman_cor <- cor(prob[,1], prob[,2], method = "spearman")

    message("bivrbetabinom: requested sigma ", sigma, ", generated ", generated_norm_cor,
            ", correlation of prob pearson ", prob_pearson_cor, ", spearman ", prob_spearman_cor)

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
if (length(args) < 6) {
  stop("Missing arguments!\nUsage: Rscript script.R <rho> <n_cells> <alpha> <beta> <mu> <size_negbinom>")
}

rho           <- as.numeric(args[1])
n             <- as.integer(args[2])
alpha_val     <- as.numeric(args[3])
beta_val      <- as.numeric(args[4])
mu_val        <- as.numeric(args[5])
size_nb_str   <- args[6] # "inf", "fixed", or numeric

config_path <- file.path("config", "config.yml")
config <- yaml::yaml.load_file(config_path)

resultsdir <- config$paths$resultsdir

alpha_A1 <- alpha_val
alpha_A2 <- alpha_val
beta_A1  <- beta_val
beta_A2  <- beta_val

set.seed(42 + n)

if (size_nb_str == "fixed") {
  sizes <- rep(mu_val, n)
} else if (size_nb_str == "inf") {
  sizes <- rpois(n, lambda = mu_val)
} else {
  size_nb_val <- as.numeric(size_nb_str)
  if (is.na(size_nb_val)) stop("Invalid size_negbinom parameter: ", size_nb_str)
  
  sizes <- rnbinom(n, size = size_nb_val, mu = mu_val)
}

sizes <- pmax(sizes, 1)

rho_str <- format(rho, nsmall = 1)
output_dir <- file.path(resultsdir, "synthetic_data", paste0("rho_", rho_str))
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
}

filename <- paste0("synthetic_counts_n", n, "_mu", mu_val, "_sizeNB", size_nb_str, "_alpha", alpha_val, "_beta", beta_val, ".tsv.gz")
output_file <- file.path(output_dir, filename)

sim <- bivrbetabinom(
  n = n,
  shape1 = c(alpha_A1, alpha_A2),
  shape2 = c(beta_A1, beta_A2),
  sizes = sizes,
  sigma = rho
)

counts <- sim$counts
x_A1_syn <- counts[,1]
x_A2_syn <- counts[,2]

sigma_used             <- sim$sigma
generated_norm_cor_val <- sim$generated_norm_cor
prob_pearson_cor_val   <- sim$prob_pearson_cor
prob_spearman_cor_val  <- sim$prob_spearman_cor

if (any(x_A1_syn > sizes) || any(x_A2_syn > sizes)) {
  stop("Generated counts violate constraint: x_A1 or x_A2 > total_reads")
}

barcodes <- paste0("CELL_", sprintf("%06d", 1:n))

out <- data.table(
  barcode = barcodes,
  total_reads = sizes,
  x_A1 = x_A1_syn,
  x_A2 = x_A2_syn
)

out[, sigma := NA_real_]
out[, generated_norm_cor := NA_real_]
out[, prob_pearson_cor := NA_real_]
out[, prob_spearman_cor := NA_real_]

out[1, `:=`(
  sigma = sigma_used,
  generated_norm_cor = generated_norm_cor_val,
  prob_pearson_cor = prob_pearson_cor_val,
  prob_spearman_cor = prob_spearman_cor_val
)]

fwrite(out, output_file, sep = "\t")
cat("Synthetic data written to:", output_file, "\n")
