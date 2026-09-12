library(data.table)
library(ggplot2)
library(patchwork)
library(ggtext)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("Usage: Rscript script.R <n> <size>")
}

CURRENT_N    <- args[1]  # e.g. "5000"
CURRENT_SIZE <- args[2]  # e.g. "1.0"

cat("Uruchamianie dla N =", CURRENT_N, "| SIZE =", CURRENT_SIZE, "\n")

config_path <- file.path("config", "config.yml")
config <- yaml::yaml.load_file(config_path)

resultsdir <- config$paths$resultsdir

INPUT_DIR       <- file.path(resultsdir, "MGLMfit_GDM_cor", "synthetic_data_extended_plus")
OUTPUT_BASE_DIR <- file.path(resultsdir, "plots_for_a_thesis", "synthetic_correlations")
INIT_VAL        <- "1e-6"
SELECTED_SEEDS  <- 0:999

RHOS <- c(-0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6)
MUS  <- c(1000, 2000, 3000, 4000, 5000)

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

format_n <- function(x) format(x, big.mark = ",", scientific = FALSE)

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

if (length(all_data_list) == 0) stop("No dtaa read. Check paths.")
dt <- rbindlist(all_data_list)


dt <- dt[as.character(n) == CURRENT_N & 
         as.character(size_nb) == CURRENT_SIZE & 
         mu %in% MUS]

dir.create(OUTPUT_BASE_DIR, recursive = TRUE, showWarnings = FALSE)


stats_dt <- dt[, .(
  n_datasets         = .N,
  mean_spearman_rho  = mean(spearman_rho, na.rm = TRUE),
  sd_spearman_rho    = sd(spearman_rho, na.rm = TRUE),
  diff_from_true_rho = mean(spearman_rho, na.rm = TRUE) - rho[1]
), by = .(rho, mu)]

setorder(stats_dt, rho, mu)

output_stats_path <- file.path(
  OUTPUT_BASE_DIR, 
  paste0("gdm_matrix_rho_vs_mu_N_", CURRENT_N, "_size_", CURRENT_SIZE, "_stats.tsv.gz")
)

fwrite(stats_dt, output_stats_path, sep = "\t", compress = "gzip")

stats_long <- melt(
  stats_dt,
  id.vars = c("rho", "mu"),
  measure.vars = c("mean_spearman_rho", "sd_spearman_rho", "diff_from_true_rho"),
  variable.name = "metric",
  value.name = "value"
)

stats_long[, rho_label := factor(
  paste0("rho = ", format(rho, nsmall = 1)),
  levels = paste0("rho = ", format(RHOS, nsmall = 1))
)]

true_rho_df <- unique(stats_long[, .(rho, rho_label)])

p_lines <- ggplot(stats_long, aes(x = mu, y = value, color = metric, group = metric)) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "gray70", linewidth = 0.4) +
  geom_hline(data = true_rho_df, aes(yintercept = rho), linetype = "dashed", color = "#555555", linewidth = 0.7) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  facet_wrap(~ rho_label, ncol = 4, scales = "fixed") +
  scale_x_continuous(breaks = MUS, labels = format_n(MUS)) +
  scale_color_manual(
    values = c(
      "mean_spearman_rho"  = "#1f77b4",
      "sd_spearman_rho"    = "#d62728",
      "diff_from_true_rho" = "#2ca02c"
    ),
    labels = c(
      "mean_spearman_rho"  = "Mean estimated Spearman's rho",
      "sd_spearman_rho"    = "SD (Standard deviation of estimated rho)",
      "diff_from_true_rho" = "Distance from a priori rho (Mean estimated rho - True a priori rho)"
    ),
    name = "Metric:"
  ) +
  labs(
    title = bquote(bold("Evaluation of estimated Spearman's correlation (N = " * .(CURRENT_N) * ", size = " * .(CURRENT_SIZE) * ")")),
    subtitle = paste0(
      "<b>mu:</b> Mean parameter of the total count per observation distribution (<i>a priori</i>) | ",
      "<b>rho:</b> <i>A priori</i> correlation used in dataset generation<br>",
      "<span style='color: #555555; font-style: italic;'>Dark dashed horizontal line indicates the true a priori rho value</span>"
    ),
    x = "mu - mean parameter of total count per observation distribution (a priori)",
    y = "Metric value"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(size = 13, hjust = 0.5, color = "#2c3e50", margin = margin(t = 6, b = 2)),
    plot.subtitle = element_markdown(size = 9.5, hjust = 0.5, color = "#2c3e50", lineheight = 1.25, margin = margin(b = 6)),
    strip.text = element_text(size = 10, face = "bold"),
    legend.position = "bottom",
    legend.title = element_text(face = "bold", size = 9.5),
    legend.text = element_text(size = 8.5),
    panel.border = element_rect(color = "#dcdde1", fill = NA, linewidth = 0.4),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8, face = "plain")
  )

output_lines_pdf_path <- file.path(
  OUTPUT_BASE_DIR, 
  paste0("gdm_lines_rho_vs_mu_N_", CURRENT_N, "_size_", CURRENT_SIZE, ".pdf")
)

pdf(output_lines_pdf_path, width = 11, height = 7)
print(p_lines)
dev.off()

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

matrix_theme <- function(row_idx, col_idx, total_rows, total_cols) {
  theme_minimal(base_size = 11) + 
    theme(
      plot.title = element_text(size = 10, face = "bold", hjust = 0.5, margin = margin(b = 2)),
      axis.title = element_blank(),
      axis.text.x = if (row_idx == total_rows) element_text(size = 9, face = "plain", color = "black") else element_blank(),
      axis.text.y = if (col_idx == 1) element_text(size = 8) else element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "#f1f2f6"),
      plot.margin = margin(t = 2, r = 5, b = 2, l = 5), 
      panel.border = element_rect(color = "#dcdde1", fill = NA, linewidth = 0.4),
      legend.position = "none"
    )
}

gdm_plots <- list()

for (r in seq_along(RHOS)) {
  for (c in seq_along(MUS)) {
    current_rho <- RHOS[r]
    current_mu  <- MUS[c]
    
    sub_dt <- dt[rho == current_rho & mu == current_mu]
    
    title_text <- if (r == 1) paste0("mu = ", format_n(current_mu)) else ""
    
    p <- ggplot(sub_dt, aes(x = spearman_rho, fill = color_group)) +
      geom_histogram(binwidth = 0.04, boundary = 0, color = "white", linewidth = 0.05, position = "stack") +
      scale_fill_manual(values = color_map, labels = legend_labels, drop = FALSE, name = "Estimation Quality:") +
      xlim(-1, 1) +
      geom_vline(xintercept = current_rho, linetype = "dashed", color = "red", linewidth = 0.6) +
      labs(title = title_text) +
      scale_y_continuous(n.breaks = 3) + 
      matrix_theme(r, c, length(RHOS), length(MUS))
    
    if (c == 1) {
      p <- p + labs(y = paste0("rho = ", format(current_rho, nsmall = 1))) + 
        theme(axis.title.y = element_text(size = 10, face = "bold", color = "#2c3e50", angle = 90, vjust = 0.5))
    }
    
    gdm_plots[[length(gdm_plots) + 1]] <- p
  }
}

combined_gdm <- wrap_plots(gdm_plots, ncol = length(MUS), nrow = length(RHOS)) +
  plot_annotation(
    title = bquote(bold("Spearman's correlation between chromatin accessibility modeled by GDM distribution (N = " * .(CURRENT_N) * ", size = " * .(CURRENT_SIZE) * ")")),
    subtitle = paste0(
      "<b>Columns:</b> <i>mu</i> - mean parameter of the total count per observation distribution (<i>a priori</i>)<br>",
      "<b>Rows:</b> <i>a priori</i> correlation (<i>rho</i>) used in dataset generation<br>",
      "<b>X-axis:</b> Estimated Spearman's correlation (recovered after modeling)<br>",
      "<b>Y-axis:</b> Number of chromatin loops per correlation bin<br>",
      "<span style='color: red; font-style: italic;'>Red dashed line indicates the true a priori rho value</span>"
    ),
    theme = theme(
      plot.title.position = "plot",
      plot.title = element_text(size = 14, hjust = 0.5, color = "#2c3e50", margin = margin(t = 6, b = 4)),
      plot.subtitle = element_markdown(size = 9.5, hjust = 0.5, color = "black", lineheight = 1.25, margin = margin(b = 6))
    )
  )

output_pdf_path <- file.path(
  OUTPUT_BASE_DIR, 
  paste0("gdm_matrix_rho_vs_mu_N_", CURRENT_N, "_size_", CURRENT_SIZE, ".pdf")
)

pdf(output_pdf_path, width = length(MUS) * 2, height = length(RHOS) + 1.5)
print(combined_gdm)
dev.off()
