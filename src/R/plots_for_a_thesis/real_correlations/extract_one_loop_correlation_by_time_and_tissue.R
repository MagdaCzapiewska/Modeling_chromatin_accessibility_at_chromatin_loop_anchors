#!/usr/bin/env Rscript

# Required libraries
required_packages <- c("data.table", "yaml", "ggplot2", "pals")
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(paste("Error: Package", pkg, "is not installed. Please install it using: install.packages('", pkg, "')"))
  }
}

library(data.table)
library(yaml)
library(ggplot2)
library(pals)

# 1. Command line arguments handling
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Error: Loop ID not provided!\nUsage: Rscript extract_loop_correlation_by_time_and_tissue.R <loop_id>")
}

target_loop <- args[1]

# 2. Load configuration and define paths
config_path <- file.path("config", "config.yml")
if (!file.exists(config_path)) {
  stop("Error: Configuration file not found at config/config.yml")
}

config <- yaml::yaml.load_file(config_path)
results_dir <- config$paths$resultsdir
pops_dir <- file.path(results_dir, "MGLMfit_GDM_cor", "real_data", "pops")

# Output folder setup
output_folder <- file.path(results_dir, "plots_for_a_thesis", "real_correlations")
dir.create(output_folder, showWarnings = FALSE, recursive = TRUE)

# 3. Dynamic search for cor_*.tsv.gz files
all_files <- list.files(pops_dir, pattern = "^cor_.*\\.tsv\\.gz$", recursive = TRUE, full.names = TRUE)
if (length(all_files) == 0) {
  cat("No result files found (cor_*.tsv.gz).\n")
  quit(status = 0)
}

res_list <- list()

# 4. Iterate over files and read data for init 1e-6
for (f in all_files) {
  init_dir <- basename(dirname(f))
  init_val <- sub("^init_", "", init_dir)
  
  if (init_val != "1e-6") next
  
  dt <- tryCatch(fread(f), error = function(e) NULL)
  if (is.null(dt) || nrow(dt) == 0) next
  if (!"loop_id" %in% names(dt)) next
  
  dt_filtered <- dt[loop_id == target_loop]
  if (nrow(dt_filtered) > 0) {
    dt_filtered[, init := init_val]
    res_list[[length(res_list) + 1]] <- dt_filtered
  }
}

# 5. Data processing
if (length(res_list) > 0) {
  final_dt <- rbindlist(res_list, fill = TRUE)
  
  # Extract tissue name
  final_dt[, tissue := sub("^[0-9]+_(.*)", "\\1", population)]
  final_dt[is.na(tissue) | tissue == "" | population == "NA" | is.na(population), tissue := "unknown"]
  
  # Define full order of time windows
  tw_early <- c("00-02", "02-04", "04-06", "06-08", "08-10")
  tw_late  <- c("10-12", "12-14", "14-16", "16-18", "18-20")
  all_tws  <- c(tw_early, tw_late)
  final_dt[, time_window := factor(time_window, levels = all_tws)]
  
  # Save original summary data as TSV to output_folder
  desired_cols <- c("time_window", "population", "tissue", "spearman_rho", "min_est_over_se")
  cols_to_keep <- intersect(desired_cols, names(final_dt))
  final_dt_print <- final_dt[, ..cols_to_keep]
  setorder(final_dt_print, time_window, population)
  fwrite(final_dt_print, file.path(output_folder, paste0("summary_", target_loop, "_init_1e-6_correlations.tsv")), sep = "\t")
  
  # Aggregation (calculating the mean correlation per tissue and time window)
  aggregated_dt <- final_dt[!is.na(spearman_rho), 
                            .(mean_spearman_rho = mean(spearman_rho, na.rm = TRUE)), 
                            by = .(time_window, tissue)]
  
  # =========================================================================
  # GLOBAL COLOR DICTIONARY (No duplicates)
  # =========================================================================
  all_unique_tissues <- sort(unique(aggregated_dt$tissue))
  alphabet_colors <- as.vector(pals::alphabet(26))
  
  if (length(all_unique_tissues) > length(alphabet_colors)) {
    alphabet_colors <- rep(alphabet_colors, length.out = length(all_unique_tissues))
  }
  strict_color_dictionary <- setNames(alphabet_colors[1:length(all_unique_tissues)], all_unique_tissues)
  # =========================================================================
  
  # --- PLOT GENERATION FUNCTION ---
  cat("\nGenerating Spearman correlation line plots (Vector PDF Output)...\n")
  
  plot_correlation_lines <- function(data_subset, title_suffix, filename) {
    if (nrow(data_subset) == 0) return(NULL)
    
    # 1. Identify unique tissues in this subset and split them in half
    active_tissues <- sort(unique(data_subset$tissue))
    total_active <- length(active_tissues)
    half_point <- ceiling(total_active / 2)
    
    # 2. Shape dictionary: 21 = Circles, 24 = Triangles
    shapes_pool <- c(rep(21, half_point), rep(24, total_active - half_point))
    shapes_dictionary <- setNames(shapes_pool, active_tissues)
    
    # Force factor levels for legend synchronization
    data_subset_plot <- copy(data_subset)
    data_subset_plot[, tissue := factor(tissue, levels = active_tissues)]
    
    # Build plot
    p <- ggplot(data_subset_plot, aes(x = time_window, y = mean_spearman_rho, group = tissue)) +
      # LINES: Colored by tissue via aes(color)
      geom_line(aes(color = tissue), size = 1.1, alpha = 0.75) +
      
      # POINTS: 
      # - color = "black" on the outside guarantees solid vector borders
      # - fill inside aes() handles the inner shape color
      # - shape inside aes() handles circle vs triangle selection
      geom_point(aes(fill = tissue, shape = tissue), color = "black", size = 3.3, stroke = 1.0, alpha = 1) + 
      
      scale_y_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.2)) +
      
      # Map dictionaries
      scale_color_manual(values = strict_color_dictionary) +
      scale_fill_manual(values = strict_color_dictionary) +
      scale_shape_manual(values = shapes_dictionary) +
      
      # Styling
      theme_minimal(base_size = 13) +
      labs(
        title = paste("Mean Spearman correlation for each tissue for loop ", target_loop, "-", title_suffix),
        x = "Time window",
        y = "Mean Spearman correlation (rho)",
        color = "Tissue type",
        fill = "Tissue type",
        shape = "Tissue type"
      ) +
      guides(
        color = guide_legend(ncol = 3, byrow = TRUE),
        fill = guide_legend(ncol = 3, byrow = TRUE),
        shape = guide_legend(ncol = 3, byrow = TRUE)
      ) +
      theme(
        plot.title = element_text(face = "bold", size = 13),
        axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
        legend.position = "right",
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_line(color = "gray92"),
        panel.grid.major.y = element_line(color = "gray92")
      )
    
    full_plot_path <- file.path(output_folder, filename)
    ggsave(full_plot_path, plot = p, width = 13, height = 5.5)
    cat("[INFO] Saved PDF plot:", filename, " (Circles:", half_point, "| Triangles:", total_active - half_point, ")\n")
  }
  
  # Plot 1: Early windows (00-10) -> Saved as PDF
  plot_correlation_lines(
    data_subset = aggregated_dt[time_window %in% tw_early], 
    title_suffix = "Windows 00-10", 
    filename = paste0("plot_", target_loop, "_spearman_lines_00_10.pdf")
  )
  
  # Plot 2: Late windows (10-20) -> Saved as PDF
  plot_correlation_lines(
    data_subset = aggregated_dt[time_window %in% tw_late], 
    title_suffix = "Windows 10-20", 
    filename = paste0("plot_", target_loop, "_spearman_lines_10_20.pdf")
  )
  
} else {
  cat(paste0("\n[!] No data found for loop ", target_loop, " with init: 1e-6.\n"))
}
