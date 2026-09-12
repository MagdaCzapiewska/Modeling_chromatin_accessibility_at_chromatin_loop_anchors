#!/usr/bin/env Rscript

library(data.table)
library(ggtern)
library(ggplot2)
library(viridis)
library(yaml)

##########################################################
# 1. PARSOWANIE ARGUMENTÓW WEJŚCIOWYCH
##########################################################

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 4) {
  stop("Użycie: Rscript plot_thesis_triangles.R <loop_id> <mode: fit|reg> <log: linear|log2> <scale_mode: global|local>")
}

loop_id_arg <- as.character(args[1])
mode_type   <- tolower(args[2])
log_arg     <- tolower(args[3])
scale_mode  <- tolower(args[4])

if (!mode_type %in% c("fit", "reg")) {
  stop("Parametr 'mode' musi przyjmować wartość 'fit' lub 'reg'.")
}

use_log2 <- log_arg %in% c("log2", "log", "true", "t", "1")
scale_label <- if (use_log2) "log2" else "linear"

if (!scale_mode %in% c("global", "local")) {
  stop("Parametr 'scale_mode' musi przyjmować wartość 'global' lub 'local'.")
}

cat(sprintf("Generowanie pojedynczych plików PDF: Loop = %s | Tryb = %s | Skala = %s | Tryb skali = %s\n",
            loop_id_arg, mode_type, scale_label, scale_mode))

##########################################################
# 2. KONFIGURACJA I ŚCIEŻKI
##########################################################

config_path <- file.path("config", "config.yml")
if (!file.exists(config_path)) {
  stop("Nie znaleziono pliku konfiguracji: ", config_path)
}
config <- yaml::yaml.load_file(config_path)
resultsdir <- config$paths$resultsdir

folder_name <- paste0("loop_", loop_id_arg, "_", mode_type, "_", scale_mode, "_", scale_label)
out_dir     <- file.path(resultsdir, "triangle_plots", "for_thesis", "time", "one_mode", folder_name)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

eval_dir <- file.path(resultsdir, "MGLMreg", "eval_0h+_time_default")

time_windows <- c("00-02", "02-04", "04-06", "06-08", "08-10", 
                  "10-12", "12-14", "14-16", "16-18", "18-20")

prefix <- paste0(mode_type, "_")
grid_res <- 150

##########################################################
# 3. FUNKCJA ANALITYCZNEJ GĘSTOŚCI GDM3
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
# 4. KROK 1: PRZYGOTOWANIE DANYCH I WYZNACZENIE SKALI GLOBALNEJ
##########################################################

plot_data_list <- list()
global_min_density <- Inf
global_max_density <- -Inf

density_legend_name <- if (use_log2) "log2(Density)" else "Density"

for (tw in time_windows) {
  eval_file <- file.path(eval_dir, paste0("eval_cor_", tw, ".tsv.gz"))
  
  if (!file.exists(eval_file)) {
    cat(sprintf("Pominięto okno %s (brak pliku eval)\n", tw))
    next
  }
  
  dt_eval <- fread(eval_file)
  dt_eval[, loop_id := as.character(loop_id)]
  row_dt <- dt_eval[loop_id == loop_id_arg]
  
  if (nrow(row_dt) == 0) {
    cat(sprintf("Brak danych dla loop_id=%s w oknie %s\n", loop_id_arg, tw))
    next
  }
  
  row_dt <- row_dt[1]
  
  a_out <- row_dt[[paste0(prefix, "alpha_x_out_est")]]
  b_out <- row_dt[[paste0(prefix, "beta_x_out_est")]]
  
  a_A1  <- row_dt[[paste0(prefix, "alpha_x_A1_est")]]
  b_A1  <- row_dt[[paste0(prefix, "beta_x_A1_est")]]
  
  a_A2  <- row_dt[[paste0(prefix, "alpha_x_A2_est")]]
  b_A2  <- row_dt[[paste0(prefix, "beta_x_A2_est")]]
  
  if (is.na(a_out) || is.na(b_out) || a_out <= 0 || b_out <= 0) {
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
  
  title_text <- sprintf("Loop: %s\nTime window: %s\nSpearman's rho = %s\n%s", 
                        loop_id_arg, tw, rho_str, metric_line)

  plot_data_list[[tw]] <- list(
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
  stop("Brak danych do wygenerowania wykresów.")
}

##########################################################
# 5. KROK 2: GENEROWANIE I ZAPIS OSOBNYCH PLIKÓW PDF
##########################################################

page_idx <- 1

for (tw in names(plot_data_list)) {
  item <- plot_data_list[[tw]]
  
  color_scale <- if (scale_mode == "global") {
    scale_color_viridis_c(
      option = "D", 
      name = density_legend_name, 
      limits = c(global_min_density, global_max_density),
      n.breaks = 5,
      guide = guide_colorbar(
        barwidth = unit(10, "cm"),
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
      # Zmniejszone paddingi wokół trójkąta i osi:
      tern.panel.expand = 0.30,
      tern.axis.arrow.sep = 0.20,
      tern.axis.text = element_text(size = 15, color = "black"),
      tern.axis.title = element_text(size = 20, face = "bold", color = "black"),
      legend.position = "bottom",
      legend.title = element_text(size = 16, face = "bold"),
      legend.margin = margin(t = -10),
      # Dociągnięcie tytułu do trójkąta poprzez ujemny dolny margines:
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
  
  pdf(single_pdf_file, width = 8, height = 8.5)
  print(p)
  dev.off()
  
  page_idx <- page_idx + 1
}

cat(sprintf("Zapisano %d indywidualnych plików PDF z ciasnym układem w folderze:\n%s\n", length(plot_data_list), out_dir))
