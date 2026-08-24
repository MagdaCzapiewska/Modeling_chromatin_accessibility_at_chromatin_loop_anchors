library(data.table)
library(ggplot2)
library(patchwork)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage: Rscript plot_synthetic_extended_gdm.R <input_dir> <expected_rho> <init_val> <output_base_dir>")
}

input_dir       <- args[1]
expected_rho    <- as.numeric(args[2])
init_val        <- args[3]
output_base_dir <- args[4]

rho_str <- format(expected_rho, nsmall = 1)

pattern <- paste0("^cor_seed.*\\.tsv\\.gz$")
files <- list.files(file.path(input_dir, paste0("rho_", rho_str), paste0("init_", init_val)), pattern = pattern, full.names = TRUE)

if (length(files) == 0) stop("No correlation files found for current rho/init.")
dt <- rbindlist(lapply(files, fread))

color_map <- c("10" = "#99FF99", "5" = "#FFFF99", "2" = "#FFCC99", "1" = "#FFB3B3", "0" = "#99CCFF", "NA" = "#D3D3D3")

param_cols <- c("alpha_x_A1_est", "alpha_x_A2_est", "alpha_x_out_est", "beta_x_A1_est", "beta_x_A2_est", "beta_x_out_est", "alpha_x_A1_SE", "alpha_x_A2_SE", "alpha_x_out_SE", "beta_x_A1_SE", "beta_x_A2_SE", "beta_x_out_SE")
dt[, n_valid_params := rowSums(!is.na(.SD)), .SDcols = param_cols]
dt[, color_group := cut(min_est_over_se, breaks = c(-Inf, 1, 2, 5, 10, Inf), labels = c("0", "1", "2", "5", "10"), right = FALSE)]
dt[n_valid_params < 8 | is.na(min_est_over_se), color_group := "NA"]
dt[, color_group := factor(color_group, levels = c("NA", "0", "1", "2", "5", "10"))]

N_CELLS <- c("500", "1000", "1500", "2000", "2500", "3000", "3500", "4000", "5000", "10000", "20000", "40000", "60000", "80000", "100000", "120000", "140000", "160000")
MUS <- c(1000, 2000, 3000, 4000, 5000, 6000)
SIZES <- c("0.1", "0.2", "0.3", "0.4", "0.5", "0.6", "0.7", "0.8", "0.9", "1.0", "1.5", "2.0", "2.5", "3.0", "3.5", "10.0", "inf", "fixed")

#N_CELLS <- c("500", "4000")
#MUS <- c(2000, 6000)
#SIZES <- c("inf", "fixed", "0.5", "10.0")

dir.create(output_base_dir, recursive = TRUE, showWarnings = FALSE)

# Helper do rysowania histogramu GDM
plot_hist_gdm <- function(sub_dt, title_txt) {
  n_datasets <- nrow(sub_dt[!is.na(spearman_rho)])
  ggplot(sub_dt, aes(x = spearman_rho, fill = color_group)) +
    geom_histogram(binwidth = 0.04, boundary = 0, color = "white", linewidth = 0.05) +
    scale_fill_manual(values = color_map, drop = FALSE, name = "Stability") +
    xlim(-1, 1) +
    geom_vline(xintercept = expected_rho, linetype = "dashed", color = "red", linewidth = 0.5) +
    labs(title = title_txt, subtitle = paste0("ds=", n_datasets), x = "Rho", y = "Count") +
    theme_minimal() + theme(plot.title = element_text(size = 7, face = "bold"), plot.subtitle = element_text(size = 6), axis.title = element_text(size = 6), legend.position = "none")
}

# --- PDF 1: Zmienia się N (Strony dla kombinacji MU x SIZE) ---
pdf(file.path(output_base_dir, "variable_N.pdf"), width = 20, height = 12)
for (current_mu in MUS) {
  for (current_size in SIZES) {
    plots <- list()
    for (current_n in N_CELLS) {
      sub_dt <- dt[mu == current_mu & size_nb == current_size & n == current_n]
      plots[[as.character(current_n)]] <- plot_hist_gdm(sub_dt, paste0("n=", current_n))
    }
    combined <- wrap_plots(plots, ncol = 6, nrow = 3) + 
      plot_annotation(title = paste0("GDM correlations | Variable N (mu=", current_mu, ", sizeNB=", current_size, ") | True rho = ", expected_rho))
    print(combined)
  }
}
dev.off()

# --- PDF 2: Zmienia się MU (Strony dla kombinacji N x SIZE) ---
pdf(file.path(output_base_dir, "variable_MU.pdf"), width = 16, height = 10)
for (current_n in N_CELLS) {
  for (current_size in SIZES) {
    plots <- list()
    for (current_mu in MUS) {
      sub_dt <- dt[n == current_n & size_nb == current_size & mu == current_mu]
      plots[[as.character(current_mu)]] <- plot_hist_gdm(sub_dt, paste0("mu=", current_mu))
    }
    combined <- wrap_plots(plots, ncol = 3, nrow = 2) + 
      plot_annotation(title = paste0("GDM correlations | Variable MU (n=", current_n, ", sizeNB=", current_size, ") | True rho = ", expected_rho))
    print(combined)
  }
}
dev.off()

# --- PDF 3: Zmienia się SIZE_NEGBINOM (Strony dla kombinacji N x MU) ---
pdf(file.path(output_base_dir, "variable_SIZE.pdf"), width = 20, height = 12)
for (current_n in N_CELLS) {
  for (current_mu in MUS) {
    plots <- list()
    for (current_size in SIZES) {
      sub_dt <- dt[n == current_n & mu == current_mu & size_nb == current_size]
      plots[[current_size]] <- plot_hist_gdm(sub_dt, paste0("size=", current_size))
    }
    combined <- wrap_plots(plots, ncol = 6, nrow = 3) + 
      plot_annotation(title = paste0("GDM correlations | Variable SizeNB (n=", current_n, ", mu=", current_mu, ") | True rho = ", expected_rho))
    print(combined)
  }
}
dev.off()
