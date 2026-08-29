#!/usr/bin/env Rscript

# Required libraries
required_packages <- c("data.table", "yaml", "ggplot2")
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(paste("Error: Package", pkg, "is not installed. Please install it using: install.packages('", pkg, "')"))
  }
}

library(data.table)
library(yaml)
library(ggplot2)

# 1. Command line arguments handling
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Error: Loop ID not provided!\nUsage: Rscript extract_one_loop_correlation_by_time_fitting_all.R <loop_id>")
}

target_loop <- args[1]

# 2. Load configuration and define paths
config_path <- file.path("config", "config.yml")
if (!file.exists(config_path)) {
  stop("Error: Configuration file not found at config/config.yml")
}

config <- yaml::yaml.load_file(config_path)
results_dir <- config$paths$resultsdir

# Input path for 'all' data
input_dir <- file.path(results_dir, "MGLMfit_GDM_cor", "real_data", "all", "init_1e-6")

# Output path
output_folder <- file.path(results_dir, "plots_for_a_thesis", "real_correlations")
dir.create(output_folder, showWarnings = FALSE, recursive = TRUE)

# 3. Define time windows order
time_windows <- c("00-02", "02-04", "04-06", "06-08", "08-10", "10-12", "12-14", "14-16", "16-18", "18-20")

# 4. Iterate over time windows files
res_list <- list()

for (tw in time_windows) {
  file_path <- file.path(input_dir, paste0("cor_", tw, ".tsv.gz"))
  
  if (file.exists(file_path)) {
    dt <- tryCatch(fread(file_path), error = function(e) NULL)
    
    if (!is.null(dt) && nrow(dt) > 0 && "loop_id" %in% names(dt)) {
      dt_sub <- dt[loop_id == target_loop]
      if (nrow(dt_sub) > 0) {
        dt_sub[, time_window := tw]
        res_list[[tw]] <- dt_sub
      }
    }
  } else {
    warning(paste("File does not exist:", file_path))
  }
}

# 5. Data processing and plot generation
if (length(res_list) > 0) {
  final_dt <- rbindlist(res_list, fill = TRUE)
  final_dt[, time_window := factor(time_window, levels = time_windows)]
  
  setorder(final_dt, time_window)
  
  # Save raw extracted correlation data as TSV
  summary_file <- file.path(output_folder, paste0("summary_", target_loop, "_init_1e-6_all_correlations.tsv"))
  fwrite(final_dt, summary_file, sep = "\t")
  cat("[INFO] Saved correlation summary TSV to:", summary_file, "\n")

  # Generate Line Plot
  cat("Generating Spearman correlation plot (Vector PDF Output)...\n")

  p <- ggplot(final_dt, aes(x = time_window, y = spearman_rho, group = 1)) +
    geom_line(color = "steelblue", size = 1.2, alpha = 0.8) +
    geom_point(color = "black", fill = "skyblue", shape = 21, size = 3.5, stroke = 1.0) +
    scale_y_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.2)) +
    theme_minimal(base_size = 13) +
    labs(
      title = paste("Spearman correlation across time windows for loop", target_loop, "(All cells)"),
      x = "Time window",
      y = "Spearman correlation (rho)"
    ) +
    theme(
      plot.title = element_text(face = "bold", size = 13),
      axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_line(color = "gray92"),
      panel.grid.major.y = element_line(color = "gray92")
    )

  full_plot_path <- file.path(output_folder, paste0("plot_", target_loop, "_spearman_lines_all.pdf"))
  ggsave(full_plot_path, plot = p, width = 10, height = 5)
  cat("[INFO] Saved PDF plot:", full_plot_path, "\n")

} else {
  cat(paste0("\n[!] No data found for loop ", target_loop, " in 'all' directory with init 1e-6.\n"))
}
