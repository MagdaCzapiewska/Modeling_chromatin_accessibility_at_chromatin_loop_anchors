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
  stop("Error: Color mode argument not provided!\n",
       "Usage: Rscript plot_all_loops_correlation_multipage.R <colored: TRUE|FALSE|colored|uncolored>")
}

# Parse boolean colored flag
color_arg <- tolower(args[1])
is_colored <- color_arg %in% c("true", "t", "1", "colored", "color")

cat("Running plot generation with is_colored =", is_colored, "\n")

# 2. Load configuration and define paths
config_path <- file.path("config", "config.yml")
if (!file.exists(config_path)) {
  stop("Error: Configuration file not found at config/config.yml")
}

config <- yaml::yaml.load_file(config_path)
datadir     <- config$paths$datadir
results_dir <- config$paths$resultsdir

# Input paths
loops_file <- file.path(datadir, "long_and_short_range_loops_D_mel.tsv")
input_dir  <- file.path(results_dir, "MGLMfit_GDM_cor", "real_data", "all", "init_1e-6")

# Output path
output_folder <- file.path(results_dir, "plots_for_a_thesis", "real_correlations")
dir.create(output_folder, showWarnings = FALSE, recursive = TRUE)

if (!file.exists(loops_file)) {
  stop("Loops file not found: ", loops_file)
}

# 3. Read loop IDs list and sort numerically
loops_table <- fread(loops_file)
if (!"loop_id" %in% names(loops_table)) {
  stop("Column 'loop_id' not found in loops file.")
}

loops_of_interest <- loops_table$loop_id
loops_of_interest <- loops_of_interest[order(as.integer(sub("^L", "", loops_of_interest)))]

# 4. Color map & legend labels definitions (if colored)
color_map <- c(
  "10" = "#99FF99", 
  "5"  = "#FFFF99", 
  "2"  = "#FFCC99", 
  "1"  = "#FFB3B3", 
  "0"  = "#99CCFF", 
  "NA" = "#D3D3D3"
)

legend_labels <- c(
  "10" = "Excellent (min EST/SE >= 10)",
  "5"  = "Good (min EST/SE >= 5)",
  "2"  = "Acceptable (min EST/SE >= 2)",
  "1"  = "Poor (min EST/SE >= 1)",
  "0"  = "Unstable (min EST/SE < 1)",
  "NA" = "Not Available / Failed"
)

# 5. Define time windows and load ALL data once into memory
time_windows <- c("00-02", "02-04", "04-06", "06-08", "08-10", "10-12", "12-14", "14-16", "16-18", "18-20")

cat("Loading data across all time windows into memory...\n")
res_list <- list()

for (tw in time_windows) {
  file_path <- file.path(input_dir, paste0("cor_", tw, ".tsv.gz"))
  
  if (file.exists(file_path)) {
    dt <- tryCatch(fread(file_path), error = function(e) NULL)
    
    if (!is.null(dt) && nrow(dt) > 0 && "loop_id" %in% names(dt)) {
      dt[, time_window := tw]
      res_list[[tw]] <- dt
    }
  } else {
    warning(paste("File does not exist:", file_path))
  }
}

if (length(res_list) == 0) {
  stop("No correlation data loaded from ", input_dir)
}

all_dt <- rbindlist(res_list, fill = TRUE)
all_dt[, time_window := factor(time_window, levels = time_windows)]

if (is_colored) {
  all_dt[, min_est_cat := fifelse(is.na(min_est_over_se), "NA",
                          fifelse(min_est_over_se >= 10, "10",
                          fifelse(min_est_over_se >= 5, "5",
                          fifelse(min_est_over_se >= 2, "2",
                          fifelse(min_est_over_se >= 1, "1", "0")))))]
  all_dt[, min_est_cat := factor(min_est_cat, levels = c("10", "5", "2", "1", "0", "NA"))]
}

setkey(all_dt, loop_id)

# 6. Multi-page PDF Plot Generation
pdf_filename <- if (is_colored) {
  "all_loops_spearman_lines_all_colored.pdf"
} else {
  "all_loops_spearman_lines_all_uncolored.pdf"
}

full_pdf_path <- file.path(output_folder, pdf_filename)

cat("Generating multi-page PDF plot at:", full_pdf_path, "...\n")

pdf(full_pdf_path, width = if (is_colored) 11 else 10, height = 5)

found_count <- 0

for (lid in loops_of_interest) {
  dt_sub <- all_dt[loop_id == lid]
  
  if (nrow(dt_sub) == 0) next
  
  dt_sub <- copy(dt_sub)
  setorder(dt_sub, time_window)
  found_count <- found_count + 1

  if (is_colored) {
    p <- ggplot(dt_sub, aes(x = time_window, y = spearman_rho, group = 1)) +
      geom_line(color = "grey40", size = 1.0, alpha = 0.8) +
      geom_point(aes(fill = min_est_cat), shape = 21, color = "black", size = 3.8, stroke = 0.8) +
      scale_y_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.2)) +
      scale_fill_manual(
        values = color_map,
        labels = legend_labels,
        drop = FALSE
      ) +
      theme_minimal(base_size = 13) +
      labs(
        title = paste("Spearman correlation across time windows for loop", lid, "(All cells)"),
        x = "Time window",
        y = "Spearman correlation (rho)",
        fill = "Estimation Quality\n(min EST / SE)"
      ) +
      theme(
        plot.title = element_text(face = "bold", size = 13),
        axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
        legend.position = "right",
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_line(color = "gray92"),
        panel.grid.major.y = element_line(color = "gray92")
      )
  } else {
    p <- ggplot(dt_sub, aes(x = time_window, y = spearman_rho, group = 1)) +
      geom_line(color = "steelblue", size = 1.2, alpha = 0.8) +
      geom_point(color = "black", fill = "skyblue", shape = 21, size = 3.5, stroke = 1.0) +
      scale_y_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.2)) +
      theme_minimal(base_size = 13) +
      labs(
        title = paste("Spearman correlation across time windows for loop", lid, "(All cells)"),
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
  }
  
  print(p)
}

dev.off()

cat("[SUCCESS] Successfully generated multi-page PDF for", found_count, "loops.\n")
cat("[INFO] Output saved to:", full_pdf_path, "\n")
