library(data.table)
library(ggplot2)
library(patchwork)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript plot_synthetic_correlations.R <input_tsv_gz> <expected_rho> <output_pdf>")
}

input_file   <- args[1]  
expected_rho <- as.numeric(args[2]) 
output_pdf   <- args[3]  

color_map <- c(
  "10" = "#99FF99",
  "5"  = "#FFFF99",
  "2"  = "#FFCC99",
  "1"  = "#FFB3B3",
  "0"  = "#99CCFF",
  "NA" = "#D3D3D3"
)

process_data <- function(dt) {
  param_cols <- c(
    "alpha_x_A1_est", "alpha_x_A2_est", "alpha_x_out_est",
    "beta_x_A1_est", "beta_x_A2_est", "beta_x_out_est",
    "alpha_x_A1_SE", "alpha_x_A2_SE", "alpha_x_out_SE",
    "beta_x_A1_SE", "beta_x_A2_SE", "beta_x_out_SE"
  )
  
  dt[, n_valid_params := rowSums(!is.na(.SD)), .SDcols = param_cols]
  
  dt[, color_group := cut(min_est_over_se, 
                          breaks = c(-Inf, 1, 2, 5, 10, Inf), 
                          labels = c("0", "1", "2", "5", "10"), 
                          right = FALSE)]
  
  dt[n_valid_params < 8 | is.na(min_est_over_se), color_group := "NA"]
  dt[, color_group := factor(color_group, levels = c("NA", "0", "1", "2", "5", "10"))]
  return(dt)
}


if (!file.exists(input_file)) {
  stop(paste("Input file not found:", input_file))
}

dt <- fread(input_file)
dt <- process_data(dt)

if(!dir.exists(dirname(output_pdf))) dir.create(dirname(output_pdf), recursive = TRUE)


pdf(output_pdf, width = 20, height = 12)


unique_n <- sort(unique(dt$n))
plots_n <- list()

for (current_n in unique_n) {
  sub_dt <- dt[n == current_n]
  n_datasets <- nrow(sub_dt[!is.na(spearman_rho)])
  
  p <- ggplot(sub_dt, aes(x = spearman_rho, fill = color_group)) +
    geom_histogram(binwidth = 0.02, boundary = 0, color = "white", linewidth = 0.1) +
    scale_fill_manual(values = color_map, drop = FALSE, name = "min_est_over_se") +
    xlim(-1, 1) +
    geom_vline(xintercept = expected_rho, linetype = "dashed", color = "red", linewidth = 0.6) +
    labs(
      title = paste0("n = ", current_n),
      subtitle = paste0("(datasets = ", n_datasets, ")"),
      x = "Spearman Rho", y = "Count"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 11, face = "bold"),
      plot.subtitle = element_text(size = 9),
      axis.title = element_text(size = 8)
    )
  plots_n[[as.character(current_n)]] <- p
}

plots_per_page <- 9
plot_chunks <- split(plots_n, ceiling(seq_along(plots_n) / plots_per_page))

for (i in seq_along(plot_chunks)) {
  page_plots <- plot_chunks[[i]]
  combined_n <- wrap_plots(page_plots, ncol = 3, byrow = TRUE) + 
    plot_layout(guides = "collect") + 
    plot_annotation(
      title = paste0("Synthetic Data Correlation Histograms (by N) | Expected Rho = ", expected_rho),
      subtitle = paste0("Page ", i, " of ", length(plot_chunks), " | Red line: true rho. Colored by stability."),
      theme = theme(plot.title = element_text(size = 14, face = "bold"))
    )
  print(combined_n)
}


unique_mu <- sort(unique(dt$mu))


plots_mu_all <- list()
for (current_mu in unique_mu) {
  sub_dt <- dt[mu == current_mu]
  n_datasets <- nrow(sub_dt[!is.na(spearman_rho)])
  
  p <- ggplot(sub_dt, aes(x = spearman_rho, fill = color_group)) +
    geom_histogram(binwidth = 0.02, boundary = 0, color = "white", linewidth = 0.1) +
    scale_fill_manual(values = color_map, drop = FALSE, name = "min_est_over_se") +
    xlim(-1, 1) +
    geom_vline(xintercept = expected_rho, linetype = "dashed", color = "red", linewidth = 0.6) +
    labs(
      title = paste0("mu = ", current_mu),
      subtitle = paste0("(all datasets = ", n_datasets, ")"),
      x = "Spearman Rho", y = "Count"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 11, face = "bold"),
      plot.subtitle = element_text(size = 9),
      axis.title = element_text(size = 8)
    )
  plots_mu_all[[as.character(current_mu)]] <- p
}

combined_mu_all <- wrap_plots(plots_mu_all, ncol = 3, nrow = 2) + 
  plot_layout(guides = "collect") + 
  plot_annotation(
    title = paste0("Synthetic Data Correlations (by MU - All Models) | Expected Rho = ", expected_rho),
    subtitle = "Red dashed line: true rho. Colored by model stability.",
    theme = theme(plot.title = element_text(size = 14, face = "bold"))
  )
print(combined_mu_all)


plots_mu_filtered <- list()
dt_stable <- dt[as.numeric(as.character(color_group)) >= 5 & !is.na(as.numeric(as.character(color_group)))]

for (current_mu in unique_mu) {
  sub_dt <- dt_stable[mu == current_mu]
  n_datasets <- nrow(sub_dt[!is.na(spearman_rho)])
  
  p <- ggplot(sub_dt, aes(x = spearman_rho, fill = color_group)) +
    geom_histogram(binwidth = 0.02, boundary = 0, color = "white", linewidth = 0.1) +
    scale_fill_manual(values = color_map, drop = FALSE, name = "min_est_over_se") +
    xlim(-1, 1) +
    geom_vline(xintercept = expected_rho, linetype = "dashed", color = "red", linewidth = 0.6) +
    labs(
      title = paste0("mu = ", current_mu),
      subtitle = paste0("(stable datasets = ", n_datasets, ")"),
      x = "Spearman Rho", y = "Count"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 11, face = "bold"),
      plot.subtitle = element_text(size = 9),
      axis.title = element_text(size = 8)
    )
  
  if (nrow(sub_dt) == 0) {
    p <- p + annotate("text", x = 0, y = 1, label = "No stable data (>=5)", size = 3, color = "red")
  }
  plots_mu_filtered[[as.character(current_mu)]] <- p
}

combined_mu_filtered <- wrap_plots(plots_mu_filtered, ncol = 3, nrow = 2) + 
  plot_layout(guides = "collect") + 
  plot_annotation(
    title = paste0("Synthetic Data Correlations (by MU - Filtered: min_est_over_se >= 5) | Expected Rho = ", expected_rho),
    subtitle = "Red dashed line: true rho. Keeping original colors for levels 5 and 10.",
    theme = theme(plot.title = element_text(size = 14, face = "bold"))
  )
print(combined_mu_filtered)


#unique_size <- unique(dt$size_nb)
#numeric_parts <- sort(as.numeric(unique_size[!unique_size %in% c("inf", "fixed")]))
#ordered_sizes <- c(unique_size[unique_size %in% c("inf", "fixed")], as.character(numeric_parts))

ordered_sizes <- c("inf", "fixed", "0.5", "1.0", "1.5", "2.0", "2.5", "3.0", "3.5", "10")

plots_size_all <- list()
for (current_size in ordered_sizes) {
  sub_dt <- dt[size_nb == current_size]
  n_datasets <- nrow(sub_dt[!is.na(spearman_rho)])
  
  p <- ggplot(sub_dt, aes(x = spearman_rho, fill = color_group)) +
    geom_histogram(binwidth = 0.02, boundary = 0, color = "white", linewidth = 0.1) +
    scale_fill_manual(values = color_map, drop = FALSE, name = "min_est_over_se") +
    xlim(-1, 1) +
    geom_vline(xintercept = expected_rho, linetype = "dashed", color = "red", linewidth = 0.6) +
    labs(
      title = paste0("sizeNB = ", current_size),
      subtitle = paste0("(all datasets = ", n_datasets, ")"),
      x = "Spearman Rho", y = "Count"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 11, face = "bold"),
      plot.subtitle = element_text(size = 9),
      axis.title = element_text(size = 8)
    )
  plots_size_all[[current_size]] <- p
}

combined_size_all <- wrap_plots(plots_size_all, ncol = 3, nrow = 4) + 
  plot_layout(guides = "collect") + 
  plot_annotation(
    title = paste0("Synthetic Data Correlations (by size_negbinom - All Models) | Expected Rho = ", expected_rho),
    subtitle = "Red dashed line: true rho. Colored by model stability.",
    theme = theme(plot.title = element_text(size = 14, face = "bold"))
  )
print(combined_size_all)

plots_size_filtered <- list()
for (current_size in ordered_sizes) {
  sub_dt <- dt_stable[size_nb == current_size]
  n_datasets <- nrow(sub_dt[!is.na(spearman_rho)])
  
  p <- ggplot(sub_dt, aes(x = spearman_rho, fill = color_group)) +
    geom_histogram(binwidth = 0.02, boundary = 0, color = "white", linewidth = 0.1) +
    scale_fill_manual(values = color_map, drop = FALSE, name = "min_est_over_se") +
    xlim(-1, 1) +
    geom_vline(xintercept = expected_rho, linetype = "dashed", color = "red", linewidth = 0.6) +
    labs(
      title = paste0("sizeNB = ", current_size),
      subtitle = paste0("(stable datasets = ", n_datasets, ")"),
      x = "Spearman Rho", y = "Count"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 11, face = "bold"),
      plot.subtitle = element_text(size = 9),
      axis.title = element_text(size = 8)
    )
  
  if (nrow(sub_dt) == 0) {
    p <- p + annotate("text", x = 0, y = 1, label = "No stable data (>=5)", size = 3, color = "red")
  }
  plots_size_filtered[[current_size]] <- p
}

combined_size_filtered <- wrap_plots(plots_size_filtered, ncol = 3, nrow = 4) + 
  plot_layout(guides = "collect") + 
  plot_annotation(
    title = paste0("Synthetic Data Correlations (by size_negbinom - Filtered: min_est_over_se >= 5) | Expected Rho = ", expected_rho),
    subtitle = "Red dashed line: true rho. Keeping original colors for levels 5 and 10.",
    theme = theme(plot.title = element_text(size = 14, face = "bold"))
  )
print(combined_size_filtered)

dev.off()
cat("Successfully generated comprehensive synthetic plot report (with optimized size_nb grid):", output_pdf, "\n")
