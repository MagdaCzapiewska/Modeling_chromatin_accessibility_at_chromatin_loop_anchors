library(data.table)
library(ggplot2)
library(patchwork)


INPUT_DIR   <- "./results/MGLMfit_GDM_cor/real_data/all/init_1e-6"
DATADIR     <- "./data"
OUTPUT_PDF  <- "./real_data_activity_matrix.pdf"

TW_LIST     <- c("06-08", "10-12", "14-16")


color_map <- c(
  "10" = "#99FF99",
  "5"  = "#FFFF99",
  "2"  = "#FFCC99",
  "1"  = "#FFB3B3",
  "0"  = "#99CCFF",
  "NA" = "#D3D3D3"
)
act_color   <- "#8A2BE2"
inact_color <- "#4B0082"

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


matrix_theme <- function(row_idx, col_idx, total_rows, total_cols) {
  theme_minimal(base_size = 11) + 
    theme(
      plot.title = element_text(size = 9, face = "bold", hjust = 0.5, margin = margin(b=1)),
      plot.subtitle = element_text(size = 7, hjust = 0.5, margin = margin(b=2)),
      axis.title = element_blank(),
      axis.text.x = if (row_idx == total_rows) element_text(size = 8, face = "bold") else element_blank(),
      axis.text.y = if (col_idx == 1) element_text(size = 8) else element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "#f1f2f6"),
      plot.margin = margin(t = 2, r = 4, b = 2, l = 4), 
      panel.border = element_rect(color = "#dcdde1", fill = NA, size = 0.4),
      legend.position = "none"
    )
}

loops_file <- file.path(DATADIR, "long_and_short_range_loops_D_mel.tsv")
if(!file.exists(loops_file)) stop("No file with loops!")
loops_data <- fread(loops_file)

master_plots <- list()


for (r in 1:2) {
  for (c in seq_along(TW_LIST)) {
    tw <- TW_LIST[c]
    path <- file.path(INPUT_DIR, paste0("cor_", tw, ".tsv.gz"))
    
    if(!file.exists(path)) {
      master_plots[[length(master_plots) + 1]] <- ggplot() + theme_void()
      master_plots[[length(master_plots) + 1]] <- ggplot() + theme_void()
      next
    }
    
    dt <- fread(path)
    dt <- process_data(dt)
    
    tw_label <- paste0(as.numeric(substr(tw, 1, 2)), "-", as.numeric(substr(tw, 4, 5)), "h")
    col_nb   <- paste0("Dmel_", tw_label, "_Neuroblasts")
    col_neur <- paste0("Dmel_", tw_label, "_Neurons")
    col_glia <- paste0("Dmel_", tw_label, "_Glia")
    
    loops_sub <- loops_data[, .(loop_id, 
                                loop_active_nb = get(col_nb), 
                                loop_active_neur = get(col_neur), 
                                loop_active_glia = get(col_glia))]
    loops_sub[, loop_active_anywhere := as.integer(loop_active_nb | loop_active_neur | loop_active_glia)]
    
    dt <- merge(dt, loops_sub, by = "loop_id", all.x = TRUE)
    
    if (r == 1) {
      sub_dt <- dt[loop_active_anywhere == 1]
    } else {
      sub_dt <- dt[loop_active_anywhere == 0]
    }
    
    n_loops_total <- nrow(sub_dt[!is.na(spearman_rho)])
    

    col_idx_hist <- (c - 1) * 2 + 1
    title_hist <- if (r == 1) paste0("time window: ", tw) else ""
    
    p_hist <- ggplot(sub_dt, aes(x = spearman_rho, fill = color_group)) +
      geom_histogram(binwidth = 0.05, boundary = 0, color = "white", linewidth = 0.05) +
      scale_fill_manual(values = color_map, drop = FALSE) +
      xlim(-1, 1) +
      labs(title = title_hist, subtitle = paste0("loops = ", n_loops_total)) +
      matrix_theme(r, col_idx_hist, 2, 6)
    
    if (col_idx_hist == 1) {
      row_label <- if (r == 1) "Anywhere Active" else "Nowhere Active"
      p_hist <- p_hist + labs(y = row_label) + 
        theme(axis.title.y = element_text(size = 10, face = "bold", color = "#2c3e50", angle = 90, vjust = 0.5))
    }
    

    col_idx_dens <- (c - 1) * 2 + 2
    title_dens <- if (r == 1) paste0("time window: ", tw) else ""
    fill_fld   <- if (r == 1) act_color else inact_color
    

    dt_stable <- sub_dt[as.numeric(as.character(color_group)) >= 5 & !is.na(as.numeric(as.character(color_group)))]
    n_loops_stable <- nrow(dt_stable[!is.na(spearman_rho)])
    
    p_dens <- ggplot(dt_stable, aes(x = spearman_rho)) +
      geom_density(fill = fill_fld, alpha = 0.5, color = "transparent") +
      xlim(-1, 1) +
      labs(title = title_dens, subtitle = paste0("loops = ", n_loops_stable)) +
      matrix_theme(r, col_idx_dens, 2, 6)
    
    if(nrow(dt_stable) < 2) {
      p_dens <- p_dens + annotate("text", x = 0, y = 0.5, label = "No data (>=5)", size = 2.5, fontface = "italic")
    }
    
    master_plots[[length(master_plots) + 1]] <- p_hist
    master_plots[[length(master_plots) + 1]] <- p_dens
  }
}

combined_matrix <- wrap_plots(master_plots, ncol = 6, nrow = 2) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = expression(bold("Real data Spearman correlation between chromatin accessibility modeled by GDM distribution")),
    theme = theme(
      plot.title.position = "plot",
      plot.title = element_text(size = 14, hjust = 0.5, color = "#2c3e50", margin = margin(t = 6, b = 2)),
      plot.subtitle = element_text(size = 10, fontface = "italic", hjust = 0.5, color = "#7f8c8d", margin = margin(b = 6)),
      legend.position = "bottom",
      legend.text = element_text(size = 9),
      legend.title = element_text(size = 9, face = "bold")
    )
  )


pdf(OUTPUT_PDF, width = 10, height = 3.5)
print(combined_matrix)
dev.off()

cat("Success! File:", OUTPUT_PDF, "\n")
