#!/usr/bin/env Rscript

# Required libraries
required_packages <- c("data.table", "yaml", "ggplot2", "scales")
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(paste("Error: Package", pkg, "is not installed. Please install it using: install.packages('", pkg, "')"))
  }
}

library(data.table)
library(yaml)
library(ggplot2)
library(scales)

# 1. Load configuration and define paths
config_path <- file.path("config", "config.yml")
if (!file.exists(config_path)) {
  stop("Error: Configuration file not found at config/config.yml")
}

config <- yaml::yaml.load_file(config_path)
resultsdir <- config$paths$resultsdir

total_reads_base_dir <- file.path(resultsdir, "counts", "total_reads")
negbinom_dir         <- file.path(resultsdir, "parameters", "negbinom_distribution_of_total_reads")

output_dir <- file.path(resultsdir, "plots_for_a_thesis", "parameters")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

time_windows <- c("00-02", "02-04", "04-06", "06-08", "08-10", 
                  "10-12", "12-14", "14-16", "16-18", "18-20")

# Custom theme for plots with larger, pure black text
custom_theme <- theme_minimal(base_size = 13) +
  theme(
    plot.title       = element_text(face = "bold", size = 16, color = "black", hjust = 0.5, margin = margin(b = 6)),
    plot.subtitle    = element_text(size = 12, color = "black", hjust = 0.5, margin = margin(b = 10)),
    axis.title       = element_text(face = "bold", size = 12, color = "black"),
    axis.text        = element_text(size = 10.5, color = "black"),
    panel.grid.minor = element_blank(),
    panel.border     = element_rect(color = "#cbd5e1", fill = NA, linewidth = 0.5)
  )

# ==============================================================================
# PART 1: NEGATIVE BINOMIAL DENSITY CURVES FOR 'ALL' POPULATION
# ==============================================================================
message("=== Generating Negative Binomial Density Curves ('all') ===")

x_seq <- 0:10000
density_list <- list()

for (tw in time_windows) {
  file_path <- file.path(negbinom_dir, paste0("negbinom_parameters_", tw, "_all_and_by_population.tsv.gz"))
  if (!file.exists(file_path)) {
    file_path <- file.path(negbinom_dir, paste0("negbinom_parameters_", tw, "_all_and_by_population.tsv"))
  }
  
  if (!file.exists(file_path)) {
    warning("Missing parameters file for time window: ", tw)
    next
  }

  dt <- fread(file_path)[population == "all"]
  if (nrow(dt) == 0 || is.na(dt$mu) || is.na(dt$size)) next

  mu <- dt$mu
  size <- dt$size
  p <- size / (size + mu)
  y_vals <- dnbinom(x_seq, size = size, prob = p)

  lbl <- paste0(tw, " (mu=", round(mu, 1), ", size=", round(size, 2), ")")
  
  density_list[[tw]] <- data.table(
    x = x_seq,
    y = y_vals,
    tw = tw,
    legend_lbl = lbl
  )
}

if (length(density_list) > 0) {
  df_density <- rbindlist(density_list)
  
  # Preserve time window ordering in legend
  unique_labels <- unique(df_density$legend_lbl)
  df_density[, legend_lbl := factor(legend_lbl, levels = unique_labels)]

  p_density <- ggplot(df_density, aes(x = x, y = y, color = legend_lbl)) +
    geom_line(linewidth = 1.1, alpha = 0.85) +
    scale_color_brewer(palette = "Spectral") +
    scale_x_continuous(labels = label_comma()) +
    labs(
      title = "Negative binomial distributions of total reads across time windows (all cells)",
      subtitle = "Comparison of fitted Probability Mass Function (PMF) curves for time windows",
      x = "Total reads per cell",
      y = "Probability Mass Function",
      color = "Time Window (Parameters)"
    ) +
    custom_theme +
    theme(
      legend.position = "right",
      legend.title    = element_text(face = "bold", size = 11.5, color = "black"),
      legend.text     = element_text(size = 10.5, color = "black")
    )

  densities_pdf <- file.path(output_dir, "negbinom_densities_all.pdf")
  densities_png <- file.path(output_dir, "negbinom_densities_all.png")

  ggsave(densities_pdf, plot = p_density, width = 11, height = 6.5)
  ggsave(densities_png, plot = p_density, width = 11, height = 6.5, dpi = 300)
  message("Saved density plot to: ", densities_pdf)
} else {
  warning("No density data was loaded.")
}

# ==============================================================================
# PART 2: EMPIRICAL TOTAL READS HISTOGRAMS FOR 'ALL' POPULATION
# ==============================================================================
message("=== Generating Empirical Total Reads Histograms ('all') ===")

raw_reads_list <- list()

for (tw in time_windows) {
  raw_reads_path <- file.path(total_reads_base_dir, paste0("total_reads_", tw, ".tsv.gz"))
  if (!file.exists(raw_reads_path)) {
    raw_reads_path <- file.path(total_reads_base_dir, paste0("total_reads_", tw, ".tsv"))
  }

  if (file.exists(raw_reads_path)) {
    raw_dt <- fread(raw_reads_path, select = "total_reads")
    raw_dt[, tw := tw]
    raw_reads_list[[tw]] <- raw_dt
  } else {
    warning("Missing raw total reads file for time window: ", tw)
  }
}

if (length(raw_reads_list) > 0) {
  df_raw <- rbindlist(raw_reads_list)

  # Calculate cell count (N) per time window and format label without extra spaces
  tw_counts <- df_raw[, .(n_cells = .N), by = tw]
  tw_counts[, tw_label := paste0("TW: ", tw, "\n(N = ", format(n_cells, big.mark = ",", trim = TRUE), ")")]

  # Ensure factor ordering follows the defined time_windows
  tw_counts[, tw := factor(tw, levels = time_windows)]
  tw_counts <- tw_counts[order(tw)]
  
  df_raw <- merge(df_raw, tw_counts[, .(tw, tw_label)], by = "tw")
  df_raw[, tw_label := factor(tw_label, levels = tw_counts$tw_label)]

  p_hist <- ggplot(df_raw, aes(x = total_reads)) +
    geom_histogram(bins = 80, fill = "#2563eb", color = "white", linewidth = 0.1, alpha = 0.85) +
    facet_wrap(~ tw_label, scales = "free_y", ncol = 3) +
    scale_x_continuous(labels = label_comma()) +
    scale_y_continuous(labels = label_comma()) +
    labs(
      title = "Empirical total reads histograms across time windows (all cells)",
      subtitle = "N - number of cells",
      x = "Total reads per cell",
      y = "Number of cells per total reads bin"
    ) +
    custom_theme +
    theme(
      legend.position  = "none",
      strip.background = element_rect(fill = "#f1f5f9", color = NA),
      strip.text       = element_text(face = "bold", size = 11, color = "black")
    )

  histograms_pdf <- file.path(output_dir, "total_reads_histograms_all.pdf")
  histograms_png <- file.path(output_dir, "total_reads_histograms_all.png")

  ggsave(histograms_pdf, plot = p_hist, width = 11, height = 9.5)
  ggsave(histograms_png, plot = p_hist, width = 11, height = 9.5, dpi = 300)
  message("Saved histograms plot to: ", histograms_pdf)
} else {
  warning("No raw total reads data was loaded.")
}

message("=== Done! Plots successfully generated in: ", output_dir, " ===")
