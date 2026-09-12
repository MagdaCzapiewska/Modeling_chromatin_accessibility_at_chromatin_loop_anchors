#!/usr/bin/env Rscript

library(data.table)
library(ggtern)
library(ggplot2)
library(viridis)
library(yaml)

##########################################################

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 5) {
  stop("Usage: Rscript plot_thesis_triangles_pop.R <loop_id> <time_window> <mode: fit|reg> <log: linear|log2> <scale_mode: global|local>")
}

loop_id_arg <- as.character(args[1])
tw_arg      <- as.character(args[2])
mode_type   <- tolower(args[3])
log_arg     <- tolower(args[4])
scale_mode  <- tolower(args[5])

if (!mode_type %in% c("fit", "reg")) {
  stop("Parameter 'mode' should be 'fit' or 'reg'.")
}

use_log2 <- log_arg %in% c("log2", "log", "true", "t", "1")
scale_label <- if (use_log2) "log2" else "linear"

if (!scale_mode %in% c("global", "local")) {
  stop("Parameter 'scale_mode' should be 'global' or 'local'.")
}

cat(sprintf("Generating PDFs for: Loop = %s | TW = %s | Mode = %s | Scale = %s | Scale mode = %s\n",
            loop_id_arg, tw_arg, mode_type, scale_label, scale_mode))

##########################################################

config_path <- file.path("config", "config.yml")
if (!file.exists(config_path)) {
  stop("File not found: ", config_path)
}
config <- yaml::yaml.load_file(config_path)
resultsdir <- config$paths$resultsdir

folder_name <- paste0("loop_", loop_id_arg, "_tw_", tw_arg, "_", mode_type, "_", scale_mode, "_", scale_label)
out_dir     <- file.path(resultsdir, "triangle_plots", "for_thesis", "time_tissue", "one_mode", folder_name)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

eval_file <- file.path(resultsdir, "MGLMreg", "eval_10h+_time_tissue_default", paste0("eval_cor_", tw_arg, ".tsv.gz"))

if (!file.exists(eval_file)) {
  stop("Eval file not found: ", eval_file)
}

dt_eval <- fread(eval_file)
dt_eval[, loop_id := as.character(loop_id)]
dt_sub <- dt_eval[loop_id == loop_id_arg]

if (nrow(dt_sub) == 0) {
  stop(sprintf("No data for loop_id=%s in time window %s", loop_id_arg, tw_arg))
}

setorder(dt_sub, population)

prefix <- paste0(mode_type, "_")
grid_res <- 150

##########################################################

dgdm3 <- function(x1, x2, x3, alpha, beta, log = FALSE, log_base = 2) {
  eps <- 1e-12
  x1 <- pmax(x1, eps)
  x2 <- pmax(x2, eps)
  x3 <- pmax(x3, eps)
  one_minus_x1 <- pmax(1 - x1, eps)

  log_B1 <- lbeta(alpha[1], beta[1])
  log_B2 <- lbeta(alpha[2], beta[2])

  log_pdf <- (alpha[1] - 1) * log(x1) +
             (alpha[2] - 1) * log(x2) +
             (beta[2] - 1) * log(x3) +
             (beta[1] - alpha[2] - beta[2]) * log(one_minus_x1) -
             log_B1 - log_B2

  if (log) {
    if (log_base == 2) {
      return(log_pdf / log(2))
    } else {
      return(log_pdf)
    }
  } else {
    return(exp(log_pdf))
  }
}

##########################################################

plot_data_list <- list()
global_min_density <- Inf
global_max_density <- -Inf

density_legend_name <- if (use_log2) "log2(Density)" else "Density"

for (i in seq_len(nrow(dt_sub))) {
  row_dt <- dt_sub[i]
  pop_id <- row_dt$population
  t_name <- if ("tissue" %in% names(row_dt) && !is.na(row_dt$tissue)) row_dt$tissue else sub("^[^_]+_", "", pop_id)

  a_out <- row_dt[[paste0(prefix, "alpha_x_out_est")]]
  b_out <- row_dt[[paste0(prefix, "beta_x_out_est")]]
  
  a_A1  <- row_dt[[paste0(prefix, "alpha_x_A1_est")]]
  b_A1  <- row_dt[[paste0(prefix, "beta_x_A1_est")]]
  
  a_A2  <- row_dt[[paste0(prefix, "alpha_x_A2_est")]]
  b_A2  <- row_dt[[paste0(prefix, "beta_x_A2_est")]]

  if (is.na(a_out) || is.na(b_out) || a_out <= 0 || b_out <= 0) {
    cat(sprintf("Skipped population %s (no parameters for x_out)\n", pop_id))
    next
  }

  if (!is.na(a_A1) && !is.na(b_A1) && a_A1 > 0 && b_A1 > 0) {
    cat2_name <- "x_A1"
    cat3_name <- "x_A2"
    a2 <- a_A1
    b2 <- b_A1
  } else if (!is.na(a_A2) && !is.na(b_A2) && a_A2 > 0 && b_A2 > 0) {
    cat2_name <- "x_A2"
    cat3_name <- "x_A1"
    a2 <- a_A2
    b2 <- b_A2
  } else {
    cat(sprintf("Skipped population %s (no second step parameters)\n", pop_id))
    next
  }

  alpha_vec <- c(a_out, a2)
  beta_vec  <- c(b_out, b2)

  mean_out <- a_out / (a_out + b_out)
  threshold <- 2 * mean_out - 1

  vals <- seq(0.0005, 0.9995, length.out = grid_res)
  grid_dt <- as.data.table(expand.grid(x_A1 = vals, x_out = vals))
  grid_dt[, x_A2 := 1 - x_A1 - x_out]
  grid_dt <- grid_dt[x_A2 > 0]

  grid_dt[, x_out_orig := x_out * (1 - threshold) + threshold]
  grid_dt[, x_A1_orig  := x_A1 * (1 - threshold)]
  grid_dt[, x_A2_orig  := x_A2 * (1 - threshold)]

  val_1 <- grid_dt$x_out_orig
  val_2 <- grid_dt[[paste0(cat2_name, "_orig")]]
  val_3 <- grid_dt[[paste0(cat3_name, "_orig")]]

  grid_dt[, density := dgdm3(val_1, val_2, val_3, alpha_vec, beta_vec, log = use_log2, log_base = 2)]

  min_d <- min(grid_dt$density, na.rm = TRUE)
  max_d <- max(grid_dt$density, na.rm = TRUE)
  if (min_d < global_min_density) global_min_density <- min_d
  if (max_d > global_max_density) global_max_density <- max_d

  breaks_LR <- pretty(c(0, 1 - threshold))
  breaks_LR <- breaks_LR[breaks_LR > 0 & breaks_LR <= 1 - threshold]
  labels_LR <- sprintf("%g", breaks_LR)

  breaks_T  <- pretty(c(threshold, 1))
  breaks_T  <- breaks_T[breaks_T >= threshold]
  labels_T  <- sprintf("%g", breaks_T)
  
  rho_val  <- row_dt[[paste0(prefix, "spearman_rho")]]
  rho_str  <- if (!is.null(rho_val) && !is.na(rho_val)) sprintf("%.2f", rho_val) else "NA"

  if (mode_type == "fit") {
    qual_val <- row_dt[["fit_min_est_over_se"]]
    qual_str <- if (!is.null(qual_val) && !is.na(qual_val)) sprintf("%.2f", qual_val) else "NA"
    metric_line <- sprintf("quality = %s", qual_str)
  } else {
    status_val <- if (!is.null(row_dt[["reg_status"]])) {
      row_dt[["reg_status"]]
    } else if (!is.null(row_dt[["reg_reg_status"]])) {
      row_dt[["reg_reg_status"]]
    } else {
      NA
    }
    status_str <- if (!is.null(status_val) && !is.na(status_val)) as.character(status_val) else "NA"
    metric_line <- sprintf("reg_status = %s", status_str)
  }
  
  title_text <- sprintf("Loop: %s\nTime window: %s\nPopulation: %s\nSpearman's rho = %s\n%s", 
                        loop_id_arg, tw_arg, pop_id, rho_str, metric_line)

  plot_data_list[[pop_id]] <- list(
    grid_dt = grid_dt,
    title_text = title_text,
    breaks_LR = breaks_LR,
    labels_LR = labels_LR,
    breaks_T = breaks_T,
    labels_T = labels_T,
    threshold = threshold
  )
}

if (length(plot_data_list) == 0) {
  stop("No data for plots.")
}

##########################################################

page_idx <- 1

for (pop_id in names(plot_data_list)) {
  item <- plot_data_list[[pop_id]]
  
  color_scale <- if (scale_mode == "global") {
    scale_color_viridis_c(
      option = "D", 
      name = density_legend_name, 
      limits = c(global_min_density, global_max_density),
      n.breaks = 5,
      guide = guide_colorbar(
        barwidth = unit(20, "cm"),
        barheight = unit(0.5, "cm"),
        title.position = "top",
        title.hjust = 0.5,
        label.theme = element_text(size = 14)
      )
    )
  } else {
    scale_color_viridis_c(
      option = "D", 
      name = density_legend_name,
      n.breaks = 5,
      guide = guide_colorbar(
        barwidth = unit(10, "cm"),
        barheight = unit(0.5, "cm"),
        title.position = "top",
        title.hjust = 0.5,
        label.theme = element_text(size = 14)
      )
    )
  }

  p <- ggtern(item$grid_dt, aes(x = x_A1, y = x_out, z = x_A2, color = density)) +
    geom_point(size = 1.3, stroke = 0) +
    labs(
      title = item$title_text,
      x = "p_A1",
      y = "p_out",
      z = "p_A2"
    ) +
    scale_L_continuous(breaks = item$breaks_LR / (1 - item$threshold), labels = item$labels_LR) +
    scale_T_continuous(breaks = (item$breaks_T - item$threshold) / (1 - item$threshold), labels = item$labels_T) +
    scale_R_continuous(breaks = item$breaks_LR / (1 - item$threshold), labels = item$labels_LR) +
    color_scale +
    theme_bw() +
    theme_showarrows() +
    theme(
      tern.panel.expand = 0.30,
      tern.axis.arrow.sep = 0.20,
      tern.axis.text = element_text(size = 15, color = "black"),
      tern.axis.title = element_text(size = 20, face = "bold", color = "black"),
      legend.position = "bottom",
      legend.title = element_text(size = 16, face = "bold"),
      legend.margin = margin(t = -10),
      plot.title = element_text(
        hjust = 0.5, 
        size = 18, 
        face = "bold", 
        lineheight = 1.15, 
        margin = margin(b = -15)
      ),
      plot.margin = margin(t = 10, r = 10, b = 10, l = 10)
    )
    
  page_num_str <- sprintf("%02d", page_idx)
  single_pdf_file <- file.path(out_dir, paste0(folder_name, "_", page_num_str, ".pdf"))
  
  pdf(single_pdf_file, width = 8, height = 9.0)
  print(p)
  dev.off()
  
  page_idx <- page_idx + 1
}

cat(sprintf("Saved %d individual PDFs for populations in directory:\n%s\n", length(plot_data_list), out_dir))
