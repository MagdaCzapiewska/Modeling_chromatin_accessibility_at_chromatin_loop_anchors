library(data.table)
library(ggplot2)
library(patchwork)

INPUT_DIR        <- "./results/MGLMfit_GDM_cor/synthetic_data_extended"
OUTPUT_BASE_DIR <- "." 
CURRENT_MU      <- 5000
CURRENT_SIZE    <- "1.0"
INIT_VAL        <- "1e-6"
SELECTED_SEEDS  <- 0:999

RHOS    <- c(-0.6, 0.0, 0.6)
N_CELLS <- c(1000, 2000, 5000, 10000, 20000)


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


all_data_list <- list()

for (r_val in RHOS) {
  rho_str <- format(r_val, nsmall = 1)
  search_path <- file.path(INPUT_DIR, paste0("rho_", rho_str), paste0("init_", INIT_VAL))
  
  files <- file.path(search_path, paste0("cor_seed", SELECTED_SEEDS, ".tsv.gz"))
  existing_files <- files[file.exists(files)]
  
  if (length(existing_files) > 0) {
    cat("Reading", length(existing_files), "files for rho =", rho_str, "\n")
    dt_rho <- rbindlist(lapply(existing_files, fread))
    all_data_list[[as.character(r_val)]] <- dt_rho
  } else {
    warning(paste("No files for rho =", rho_str, "on path:", search_path))
  }
}

if (length(all_data_list) == 0) stop("No data read. Check paths.")
dt <- rbindlist(all_data_list)


dt <- dt[mu == CURRENT_MU & size_nb == CURRENT_SIZE]


param_cols <- c(
  "alpha_x_A1_est", "alpha_x_A2_est", "alpha_x_out_est", 
  "beta_x_A1_est", "beta_x_A2_est", "beta_x_out_est", 
  "alpha_x_A1_SE", "alpha_x_A2_SE", "alpha_x_out_SE", 
  "beta_x_A1_SE", "beta_x_A2_SE", "beta_x_out_SE"
)

dt[, n_valid_params := rowSums(!is.na(.SD)), .SDcols = param_cols]
dt[, color_group := cut(min_est_over_se, breaks = c(-Inf, 1, 2, 5, 10, Inf), 
                        labels = c("0", "1", "2", "5", "10"), right = FALSE)]
dt[n_valid_params < 8 | is.na(min_est_over_se), color_group := "NA"]
dt[, color_group := factor(color_group, levels = c("NA", "0", "1", "2", "5", "10"))]

dir.create(OUTPUT_BASE_DIR, recursive = TRUE, showWarnings = FALSE)


matrix_theme <- function(row_idx, col_idx, total_rows, total_cols) {
  theme_minimal(base_size = 11) + 
    theme(
      plot.title = element_text(size = 10, face = "bold", hjust = 0.5, margin = margin(b=2)),
      axis.title = element_blank(),
      axis.text.x = if (row_idx == total_rows) element_text(size = 9, face = "bold", color = "black") else element_blank(),
      axis.text.y = if (col_idx == 1) element_text(size = 8) else element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "#f1f2f6"),
      plot.margin = margin(t = 2, r = 5, b = 2, l = 5), 
      panel.border = element_rect(color = "#dcdde1", fill = NA, linewidth = 0.4),
      legend.position = "none"
    )
}


cat("Generating PDF... \n")
gdm_plots <- list()

for (r in seq_along(RHOS)) {
  for (c in seq_along(N_CELLS)) {
    current_rho <- RHOS[r]
    current_n   <- as.integer(N_CELLS[c])
    
    sub_dt <- dt[rho == current_rho & n == current_n]
    
    title_text <- if (r == 1) paste0("N = ", current_n) else ""
    
    p <- ggplot(sub_dt, aes(x = spearman_rho, fill = color_group)) +
      geom_histogram(binwidth = 0.04, boundary = 0, color = "white", linewidth = 0.05, position = "stack") +
      scale_fill_manual(values = color_map, labels = legend_labels, drop = FALSE, name = "Estimation Quality:") +
      xlim(-1, 1) +
      geom_vline(xintercept = current_rho, linetype = "dashed", color = "red", linewidth = 0.6) +
      labs(title = title_text) +
      scale_y_continuous(n.breaks = 3) + 
      matrix_theme(r, c, length(RHOS), length(N_CELLS))
    
    if (c == 1) {
      p <- p + labs(y = paste0("rho = ", format(current_rho, nsmall = 1))) + 
        theme(axis.title.y = element_text(size = 10, face = "bold", color = "#2c3e50", angle = 90, vjust = 0.5))
    }
    
    gdm_plots[[length(gdm_plots) + 1]] <- p
  }
}


combined_gdm <- wrap_plots(gdm_plots, ncol = 5, nrow = 3) +
  plot_annotation(
    title = expression(bold("Spearman correlation between chromatin accessibility modeled by GDM distribution (" * mu * " = 5000, size = 1.0)")),
    subtitle = expression(italic("Red dashed line indicates the true simulated correlation value (rho = -0.6, 0.0, or 0.6)")),
    theme = theme(
      plot.title.position = "plot",
      plot.title = element_text(size = 14, hjust = 0.5, color = "#2c3e50", margin = margin(t = 6, b = 2)),
      plot.subtitle = element_text(size = 10, fontface = "italic", hjust = 0.5, color = "red", margin = margin(b = 6))
    )
  )


pdf(file.path(OUTPUT_BASE_DIR, "gdm_matrix_by_stability.pdf"), width = 10, height = 4.0)
print(combined_gdm)
dev.off()

cat("Success.\n")
