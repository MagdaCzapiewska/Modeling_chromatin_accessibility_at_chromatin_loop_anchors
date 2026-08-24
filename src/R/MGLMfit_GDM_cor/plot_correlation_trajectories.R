library(data.table)
library(ggplot2)
library(patchwork)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript plot_correlation_trajectories.R <type> <input_dir> <output_pdf>")
}

data_type  <- args[1]  # "gdm" or "naive"
input_dir  <- args[2]  # Directory containing cor_{tw}.tsv.gz files
output_pdf <- args[3]  # Output PDF path

tw_list <- c("00-02", "02-04", "04-06", "06-08", "08-10", "10-12", "12-14", "14-16", "16-18", "18-20")

# ==============================================================================
# 1. DATA LOADING & PREPARATION
# ==============================================================================
dt_list <- list()

cat("Loading all time windows...\n")
for (tw in tw_list) {
  path <- file.path(input_dir, paste0("cor_", tw, ".tsv.gz"))
  if (!file.exists(path)) {
    warning(paste("Missing file:", path))
    next
  }
  
  dt_sub <- fread(path)
  
  # Generate unique loop identifier if missing
  if (!"loop_id" %in% names(dt_sub)) {
    if ("chr_x" %in% names(dt_sub)) {
      dt_sub[, loop_id := paste(chr_x, start_x, chr_y, start_y, sep = "_")]
    } else {
      dt_sub[, loop_id := paste0("loop_", .I)]
    }
  }
  
  # Select minimal columns
  cols_to_keep <- c("loop_id", "spearman_rho", "pearson_r")
  dt_minimal <- dt_sub[, ..cols_to_keep]
  dt_minimal[, tw_window := tw]
  dt_list[[tw]] <- dt_minimal
}

main_dt <- rbindlist(dt_list)
main_dt[, tw_window := factor(tw_window, levels = tw_list)]

melted_dt <- melt(main_dt, id.vars = c("loop_id", "tw_window"), 
                  measure.vars = c("spearman_rho", "pearson_r"),
                  variable.name = "metric", value.name = "value")

melted_dt <- melted_dt[!is.na(value)]
total_loops <- uniqueN(melted_dt$loop_id)

# ==============================================================================
# 2. PLOT GENERATION (ALL LOOPS, OPAQUE LINES)
# ==============================================================================
if(!dir.exists(dirname(output_pdf))) dir.create(dirname(output_pdf), recursive = TRUE)

# Wide PDF device
pdf(output_pdf, width = 18, height = 11)

cat("Rendering raw trajectory lines for all loops...\n")
cat("Note: This will likely result in heavy overplotting.\n")

# HYPERPARAMETERS FOR OPAQUE LINES
line_alpha <- 1.0  # FULLY OPAQUE
line_width <- 0.3  # Standard line width

# --- SPEARMAN PANEL ---
p_spearman <- ggplot(melted_dt[metric == "spearman_rho"], aes(x = tw_window, y = value, group = loop_id)) +
  # Draw every single loop profile as an opaque line
  geom_line(alpha = line_alpha, color = "#4682B4", linewidth = line_width) +
  ylim(-1, 1) +
  labs(
    title = "Raw Profiles: Spearman Rho",
    subtitle = paste0("Plotting ALL unique loops (opaque lines) | n = ", total_loops),
    x = "Time Window", y = "Spearman Rho"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    panel.grid.minor = element_blank()
  )

# --- PEARSON PANEL ---
p_pearson <- ggplot(melted_dt[metric == "pearson_r"], aes(x = tw_window, y = value, group = loop_id)) +
  # Draw every single loop profile as an opaque line
  geom_line(alpha = line_alpha, color = "#708090", linewidth = line_width) +
  ylim(-1, 1) +
  labs(
    title = "Raw Profiles: Pearson R",
    subtitle = paste0("Plotting ALL unique loops (opaque lines) | n = ", total_loops),
    x = "Time Window", y = "Pearson R"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.title.y = element_blank(),
    panel.grid.minor = element_blank()
  )

# Combine the panels
final_plot <- p_spearman + p_pearson + 
  plot_annotation(
    title = paste("Raw Chromatin Loop Correlation Trajectories -", toupper(data_type), "Dataset"),
    subtitle = paste0("Visualizing the complete pathway of every individual loop without aggregation."),
    theme = theme(plot.title = element_text(size = 18, face = "bold"), plot.subtitle = element_text(size = 12))
  )

print(final_plot)

dev.off()
cat("Successfully exported raw trajectory field to:", output_pdf, "\n")
