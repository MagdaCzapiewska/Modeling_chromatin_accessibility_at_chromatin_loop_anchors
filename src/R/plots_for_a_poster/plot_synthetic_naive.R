library(data.table)
library(ggplot2)
library(patchwork)

INPUT_DIR        <- "./results/naive_correlation/synthetic_data_extended"
OUTPUT_BASE_DIR <- "." 
EXPECTED_RHO    <- 0.0
CURRENT_MU      <- 5000
SELECTED_SEEDS  <- 0:999

N_CELLS <- c(1000, 2000, 5000, 10000, 20000)
SIZES   <- c("0.5", "1.0", "2.0", "inf", "fixed")


rho_str <- format(EXPECTED_RHO, nsmall = 1)
search_path <- file.path(INPUT_DIR, paste0("rho_", rho_str))
files <- file.path(search_path, paste0("naive_cor_seed", SELECTED_SEEDS, ".tsv.gz"))
existing_files <- files[file.exists(files)]

if (length(existing_files) == 0) stop("No files with data!")
dt <- rbindlist(lapply(existing_files, fread))

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
      panel.border = element_rect(color = "#dcdde1", fill = NA, linewidth = 0.4)
    )
}

cat("Generating PDF for Pearson... \n")
pearson_plots <- list()

for (r in seq_along(SIZES)) {
  for (c in seq_along(N_CELLS)) {
    current_size <- SIZES[r]
    current_n    <- N_CELLS[c]
    
    sub_dt <- dt[mu == CURRENT_MU & size_nb == current_size & n == current_n]
    
    title_text <- if (r == 1) paste0("N = ", current_n) else ""
    
    p <- ggplot(sub_dt, aes(x = pearson_r)) +
      geom_histogram(binwidth = 0.04, boundary = 0, fill = "#4682B4", color = "white", linewidth = 0.05) +
      xlim(-0.4, 0.4) +
      geom_vline(xintercept = EXPECTED_RHO, linetype = "dashed", color = "#e74c3c", linewidth = 0.6) +
      labs(title = title_text) +
      scale_y_continuous(n.breaks = 3) + 
      matrix_theme(r, c, length(SIZES), length(N_CELLS))
    
    if (c == 1) {
      row_label <- ifelse(current_size == "fixed", "fixed", paste0("s = ", current_size))
      p <- p + labs(y = row_label) + 
        theme(axis.title.y = element_text(size = 10, face = "bold", color = "#2c3e50", angle = 90, vjust = 0.5))
    }
    pearson_plots[[length(pearson_plots) + 1]] <- p
  }
}

combined_pearson <- wrap_plots(pearson_plots, ncol = 5, nrow = 5) +
  plot_annotation(
    title = expression(bold("Naive Pearson correlation between counts at chromatin loop anchors (mu = 5000, size = s)")),
    subtitle = "Red dashed line indicates the true expected correlation (r = 0.0)",
    theme = theme(
      plot.title.position = "plot",
      plot.title = element_text(size = 14, hjust = 0.5, color = "#2c3e50", margin = margin(t = 6, b = 2)),
      plot.subtitle = element_text(size = 10, fontface = "italic", hjust = 0.5, color = "#e74c3c", margin = margin(b = 8))
    )
  )

pdf(file.path(OUTPUT_BASE_DIR, "naive_matrix_pearson.pdf"), width = 10, height = 4)
print(combined_pearson)
dev.off()

cat("Generating PDF for Spearman... \n")
spearman_plots <- list()

for (r in seq_along(SIZES)) {
  for (c in seq_along(N_CELLS)) {
    current_size <- SIZES[r]
    current_n    <- N_CELLS[c]
    
    sub_dt <- dt[mu == CURRENT_MU & size_nb == current_size & n == current_n]
    
    title_text <- if (r == 1) paste0("N = ", current_n) else ""
    
    p <- ggplot(sub_dt, aes(x = spearman_rho)) +
      geom_histogram(binwidth = 0.04, boundary = 0, fill = "#708090", color = "white", linewidth = 0.05) +
      xlim(-0.4, 0.4) +
      geom_vline(xintercept = EXPECTED_RHO, linetype = "dashed", color = "#e74c3c", linewidth = 0.6) +
      labs(title = title_text) +
      scale_y_continuous(n.breaks = 3) + 
      matrix_theme(r, c, length(SIZES), length(N_CELLS))
    
    if (c == 1) {
      row_label <- ifelse(current_size == "fixed", "fixed", paste0("s = ", current_size))
      p <- p + labs(y = row_label) + 
        theme(axis.title.y = element_text(size = 10, face = "bold", color = "#2c3e50", angle = 90, vjust = 0.5))
    }
    spearman_plots[[length(spearman_plots) + 1]] <- p
  }
}

combined_spearman <- wrap_plots(spearman_plots, ncol = 5, nrow = 5) +
  plot_annotation(
    title = expression(bold("Correlation between counts at chromatin loop anchors (mu = 5000, size = s)")),
    subtitle = expression(italic("Red dashed line indicates the true expected correlation (" * rho * " = 0.0)")),
    theme = theme(
      plot.title.position = "plot",
      plot.title = element_text(size = 14, hjust = 0.5, color = "#2c3e50", margin = margin(t = 6, b = 2)),
      plot.subtitle = element_text(size = 10, hjust = 0.5, color = "#e74c3c", margin = margin(b = 8))
    )
  )

pdf(file.path(OUTPUT_BASE_DIR, "naive_matrix_spearman.pdf"), width = 10, height = 4)
print(combined_spearman)
dev.off()

cat("Success.\n")
