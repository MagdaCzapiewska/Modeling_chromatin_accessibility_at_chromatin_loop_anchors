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

# 1. Load configuration
config_path <- file.path("config", "config.yml")
if (!file.exists(config_path)) {
  stop("Error: Configuration file not found at config/config.yml")
}

config <- yaml::yaml.load_file(config_path)
results_dir <- config$paths$resultsdir
datadir <- config$paths$datadir

loops_file <- file.path(datadir, "long_and_short_range_loops_D_mel.tsv")

if (!file.exists(loops_file)) {
  stop("Loops file not found: ", loops_file)
}

# 2. Read loop IDs list and sort numerically
loops_table <- fread(loops_file)
if (!"loop_id" %in% names(loops_table)) {
  stop("Column 'loop_id' not found in loops file.")
}

loops_of_interest <- loops_table$loop_id
loops_of_interest <- loops_of_interest[order(as.integer(sub("^L", "", loops_of_interest)))]
num_loops <- length(loops_of_interest) # 417 loops

# 3. Define paths and time windows
all_dir <- file.path(results_dir, "MGLMfit_GDM_cor", "real_data", "all", "init_1e-6")
pops_dir <- file.path(results_dir, "MGLMfit_GDM_cor", "real_data", "pops", "init_1e-6")

output_folder <- file.path(results_dir, "plots_for_a_thesis", "real_correlations", "quality")
dir.create(output_folder, showWarnings = FALSE, recursive = TRUE)

time_windows <- c("00-02", "02-04", "04-06", "06-08", "08-10", "10-12", "12-14", "14-16", "16-18", "18-20")

# 4. Color map definitions
color_map <- c(
  "10" = "#99FF99", 
  "5"  = "#FFFF99", 
  "2"  = "#FFCC99", 
  "1"  = "#FFB3B3", 
  "0"  = "#99CCFF", 
  "NA" = "#D3D3D3"
)

# Categorization helper function
categorize_quality <- function(dt) {
  dt[, min_est_cat := fifelse(is.na(min_est_over_se), "NA",
                      fifelse(min_est_over_se >= 10, "10",
                      fifelse(min_est_over_se >= 5, "5",
                      fifelse(min_est_over_se >= 2, "2",
                      fifelse(min_est_over_se >= 1, "1", "0")))))]
  
  dt[, min_est_cat := factor(min_est_cat, levels = c("NA", "0", "1", "2", "5", "10"))]
  return(dt)
}

# Helper function to generate x-axis labels with sample size (n = ...)
get_hq_x_labels <- function(dt_hq, time_windows) {
  counts_dt <- dt_hq[, .(n = .N), by = time_window]
  full_counts <- data.table(time_window = time_windows)
  full_counts <- merge(full_counts, counts_dt, by = "time_window", all.x = TRUE)
  full_counts[is.na(n), n := 0]
  
  labels_vec <- setNames(
    paste0(full_counts$time_window, "\n(n_fits=", format(full_counts$n, big.mark = ",", , trim = TRUE), ")"),
    full_counts$time_window
  )
  return(labels_vec[time_windows])
}

# Unified custom theme for barplots
custom_theme <- theme_minimal(base_size = 12) +
  theme(
    plot.title         = element_text(size = 13, face = "bold", hjust = 0.5, color = "black", margin = margin(b = 6)),
    plot.subtitle      = element_text(size = 11, hjust = 0.5, color = "black", margin = margin(b = 8), lineheight = 1.2),
    axis.title         = element_text(size = 11, face = "bold", color = "black"),
    axis.text          = element_text(size = 10, color = "black"),
    axis.text.x        = element_text(angle = 45, hjust = 1, color = "black"),
    legend.position    = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.border       = element_rect(color = "#dcdde1", fill = NA, linewidth = 0.5)
  )

# Theme for violin plots (centered x-axis labels)
violin_theme <- custom_theme + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 9.5))

# ==============================================================================
# PART 1: ALL CELLS DATA
# ==============================================================================

master_all_grid <- CJ(
  time_window = time_windows,
  loop_id = loops_of_interest
)

all_list <- list()
for (tw in time_windows) {
  fpath <- file.path(all_dir, paste0("cor_", tw, ".tsv.gz"))
  if (!file.exists(fpath)) {
    fpath <- file.path(all_dir, paste0("cor_", tw, ".tsv"))
  }
  
  if (file.exists(fpath)) {
    dt <- fread(fpath)
    dt[, time_window := tw]
    all_list[[tw]] <- dt
  } else {
    warning("File not found for 'all' time window: ", tw)
  }
}

dt_all_raw <- rbindlist(all_list, fill = TRUE)

dt_all <- merge(master_all_grid, dt_all_raw, by = c("time_window", "loop_id"), all.x = TRUE)
dt_all <- categorize_quality(dt_all)

# Statistics table for 'all'
stats_all <- dt_all[, .(count = .N), by = .(time_window, min_est_cat)]

full_stats_all_grid <- CJ(
  time_window = time_windows,
  min_est_cat = factor(c("10", "5", "2", "1", "0", "NA"), levels = c("10", "5", "2", "1", "0", "NA"))
)

stats_all <- merge(full_stats_all_grid, stats_all, by = c("time_window", "min_est_cat"), all.x = TRUE)
stats_all[is.na(count), count := 0]
stats_all[, total_loops := sum(count), by = time_window]
stats_all[, percentage := round((count / total_loops) * 100, 2)]

fwrite(stats_all, file.path(output_folder, "mglm_quality_stats_all.tsv"), sep = "\t")

# Barplot for 'all'
p_all <- ggplot(dt_all, aes(x = factor(time_window, levels = time_windows), fill = min_est_cat)) +
  geom_bar(position = "stack", width = 0.7, color = "black", linewidth = 0.2) +
  scale_fill_manual(values = color_map, drop = FALSE) +
  scale_y_continuous(labels = label_comma()) +
  labs(
    title = "Fit quality across time windows (all cells)",
    subtitle = paste0(
      "Fit for each combination of loop and time window (all cells)\n",
      "Total number of loops: ", format(num_loops, big.mark = ","), "\n",
      "Total number of fits: ", format(nrow(master_all_grid), big.mark = ",")
    ),
    x = "Time window",
    y = "Number of loops"
  ) +
  custom_theme

ggsave(file.path(output_folder, "mglm_quality_all.pdf"), plot = p_all, width = 8, height = 6)
ggsave(file.path(output_folder, "mglm_quality_all.png"), plot = p_all, width = 8, height = 6, dpi = 300)

# --- Violin plot for 'all' (min_est_over_se >= 5) ---
dt_all_hq <- dt_all[!is.na(min_est_over_se) & min_est_over_se >= 5]

if (nrow(dt_all_hq) > 0) {
  x_labels_all <- get_hq_x_labels(dt_all_hq, time_windows)
  dt_all_hq[, time_window_f := factor(time_window, levels = time_windows)]
  
  p_all_violin <- ggplot(dt_all_hq, aes(x = time_window_f, y = spearman_rho)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.5) +
    geom_violin(fill = "#74b9ff", color = "#0984e3", alpha = 0.7, linewidth = 0.4, trim = FALSE) +
    geom_boxplot(width = 0.12, fill = "white", color = "black", outlier.size = 0.4, linewidth = 0.4) +
    scale_x_discrete(labels = x_labels_all) +
    scale_y_continuous(limits = c(-1, 1), breaks = seq(-1, 1, by = 0.2)) +
    labs(
      title = "Estimated correlation distribution across time windows (fits for all cells)",
      subtitle = paste0("High-quality fits only (min_est_over_se >= 5)\nTotal high-quality fits: ", format(nrow(dt_all_hq), big.mark = ",")),
      x = "Time window",
      y = "Estimated Spearman's correlation"
    ) +
    violin_theme

  ggsave(file.path(output_folder, "mglm_cor_violin_all.pdf"), plot = p_all_violin, width = 9, height = 6)
  ggsave(file.path(output_folder, "mglm_cor_violin_all.png"), plot = p_all_violin, width = 9, height = 6, dpi = 300)
}

# ==============================================================================
# PART 2: POPULATIONS DATA
# ==============================================================================

pop_files <- list.files(pops_dir, pattern = "^cor_.*\\.tsv(\\.gz)?$", recursive = TRUE, full.names = TRUE)

pops_list <- list()
for (fpath in pop_files) {
  dt <- fread(fpath)
  
  if (!"time_window" %in% names(dt)) {
    fname <- basename(fpath)
    tw_match <- regmatches(fname, regexpr("\\d{2}-\\d{2}", fname))
    if (length(tw_match) > 0) {
      dt[, time_window := tw_match]
    }
  }
  pops_list[[fpath]] <- dt
}

dt_pops_raw <- rbindlist(pops_list, fill = TRUE)

unique_pop_units <- unique(dt_pops_raw[!is.na(population) & !is.na(time_window), .(time_window, population)])
num_pop_units <- nrow(unique_pop_units)

cat("Identified", num_pop_units, "unique population units (time_window + population pairs) across all time windows.\n")

master_pops_grid <- as.data.table(
  cbind(
    unique_pop_units[rep(seq_len(.N), each = num_loops)],
    loop_id = rep(loops_of_interest, times = num_pop_units)
  )
)

dt_pops <- merge(master_pops_grid, dt_pops_raw, by = c("loop_id", "time_window", "population"), all.x = TRUE)
dt_pops <- categorize_quality(dt_pops)

# Statistics table for 'pops'
stats_pops <- dt_pops[, .(count = .N), by = .(time_window, min_est_cat)]

full_stats_pops_grid <- CJ(
  time_window = time_windows,
  min_est_cat = factor(c("10", "5", "2", "1", "0", "NA"), levels = c("10", "5", "2", "1", "0", "NA"))
)

stats_pops <- merge(full_stats_pops_grid, stats_pops, by = c("time_window", "min_est_cat"), all.x = TRUE)
stats_pops[is.na(count), count := 0]
stats_pops[, total_fits := sum(count), by = time_window]
stats_pops[, percentage := round((count / total_fits) * 100, 2)]

fwrite(stats_pops, file.path(output_folder, "mglm_quality_stats_pops.tsv"), sep = "\t")

# Barplot for 'pops'
p_pops <- ggplot(dt_pops, aes(x = factor(time_window, levels = time_windows), fill = min_est_cat)) +
  geom_bar(position = "stack", width = 0.7, color = "black", linewidth = 0.2) +
  scale_fill_manual(values = color_map, drop = FALSE) +
  scale_y_continuous(labels = label_comma()) +
  labs(
    title = "Fit quality across time windows (populations)",
    subtitle = paste0(
      "Fit for each combination of loop and population\n",
      "Total number of loops: ", format(num_loops, big.mark = ","), "\n",
      "Total number of populations across all time windows: ", format(num_pop_units, big.mark = ","), "\n",
      "Total number of fits: ", format(nrow(master_pops_grid), big.mark = ",")
    ),
    x = "Time window",
    y = "Number of fits"
  ) +
  custom_theme

ggsave(file.path(output_folder, "mglm_quality_pops.pdf"), plot = p_pops, width = 8, height = 6.5)
ggsave(file.path(output_folder, "mglm_quality_pops.png"), plot = p_pops, width = 8, height = 6.5, dpi = 300)

# --- Violin plot for 'pops' (min_est_over_se >= 5) ---
dt_pops_hq <- dt_pops[!is.na(min_est_over_se) & min_est_over_se >= 5]

if (nrow(dt_pops_hq) > 0) {
  x_labels_pops <- get_hq_x_labels(dt_pops_hq, time_windows)
  dt_pops_hq[, time_window_f := factor(time_window, levels = time_windows)]
  
  p_pops_violin <- ggplot(dt_pops_hq, aes(x = time_window_f, y = spearman_rho)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.5) +
    geom_violin(fill = "#55efc4", color = "#00b894", alpha = 0.7, linewidth = 0.4, trim = FALSE) +
    geom_boxplot(width = 0.12, fill = "white", color = "black", outlier.size = 0.4, linewidth = 0.4) +
    scale_x_discrete(labels = x_labels_pops) +
    scale_y_continuous(limits = c(-1, 1), breaks = seq(-1, 1, by = 0.2)) +
    labs(
      title = "Estimated correlation distribution across time windows (fits for populations)",
      subtitle = paste0("High-quality fits only (min_est_over_se >= 5)\nTotal high-quality fits: ", format(nrow(dt_pops_hq), big.mark = ",")),
      x = "Time window",
      y = "Estimated Spearman's correlation"
    ) +
    violin_theme

  ggsave(file.path(output_folder, "mglm_cor_violin_pops.pdf"), plot = p_pops_violin, width = 9, height = 6.5)
  ggsave(file.path(output_folder, "mglm_cor_violin_pops.png"), plot = p_pops_violin, width = 9, height = 6.5, dpi = 300)
}

cat("Successfully generated statistics, quality barplots, and correlation violin plots in:", output_folder, "\n")
