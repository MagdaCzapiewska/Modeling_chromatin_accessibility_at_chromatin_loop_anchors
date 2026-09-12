#!/usr/bin/env Rscript

library(data.table)
library(ggtern)
library(ggplot2)
library(viridis)
library(yaml)

##########################################################
# PARSOWANIE ARGUMENTÓW WEJŚCIOWYCH
##########################################################

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("Użycie: Rscript plot_loop_triangles_dens.R <loop_id> <fit|reg> [linear|log2]")
}

loop_id_arg <- as.character(args[1])
mode_type   <- tolower(args[2]) # "fit" lub "reg"
scale_type  <- if (length(args) >= 3) tolower(args[3]) else "linear"

if (!mode_type %in% c("fit", "reg")) {
  stop("Parametr 'mode_type' musi przyjmować wartość 'fit' lub 'reg'.")
}

use_log2 <- scale_type %in% c("log2", "log", "true", "t", "1")
scale_label <- if (use_log2) "log2" else "linear"

cat(sprintf("Uruchamianie generowania gęstości: Loop = %s | Tryb = %s | Skala = %s\n", 
            loop_id_arg, mode_type, scale_label))

##########################################################
# KONFIGURACJA I ŚCIEŻKI
##########################################################

config_path <- file.path("config", "config.yml")
if (!file.exists(config_path)) {
  stop("Nie znaleziono pliku konfiguracji: ", config_path)
}
config <- yaml::yaml.load_file(config_path)
resultsdir <- config$paths$resultsdir

out_dir <- file.path(resultsdir, "triangle_plots", "time", paste0("density_", mode_type))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
pdf_file <- file.path(out_dir, paste0("loop_", loop_id_arg, "_density_combined_", scale_label, ".pdf"))

eval_dir <- file.path(resultsdir, "MGLMreg", "eval_0h+_time_default")

time_windows <- c("00-02", "02-04", "04-06", "06-08", "08-10", 
                  "10-12", "12-14", "14-16", "16-18", "18-20")

prefix <- paste0(mode_type, "_")
grid_res <- 150

##########################################################
# FUNKCJA ANALITYCZNEJ GĘSTOŚCI GDM3 (Z OBSŁUGĄ LOG2)
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
# KROK 1: OBLICZENIE GĘSTOŚCI I ZAKRESU GLOBALNEGO
##########################################################

plot_data_list <- list()
global_min_density <- Inf
global_max_density <- -Inf

density_legend_name <- if (use_log2) "$\\log_2(\\text{Density})$" else "Density"

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
  labels_LR <- paste0("$", sub("-0", "-", sub("e(.*)", " \\\\cdot 10^{\\1}", labels_LR)), "$")

  breaks_T  <- pretty(c(threshold, 1))
  breaks_T  <- breaks_T[breaks_T >= threshold]
  labels_T  <- sub(" - 0$", "", sprintf("1 - %g", 1 - breaks_T))
  labels_T  <- paste0("$", sub("-0", "-", sub("e(.*)", " \\\\cdot 10^{\\1}", labels_T)), "$")
  
  rho_val <- row_dt[[paste0(prefix, "spearman_rho")]]
  rho_str <- if (!is.null(rho_val) && !is.na(rho_val)) sprintf("%.3f", rho_val) else "NA"
  
  title_text <- sprintf("loop_id: %s | time_window: %s", loop_id_arg, tw)
  
  if (mode_type == "fit") {
    qual_val <- row_dt[["fit_min_est_over_se"]]
    qual_str <- if (!is.null(qual_val) && !is.na(qual_val)) sprintf("%.3f", qual_val) else "NA"
    subtitle_text <- sprintf("Spearman's $\\rho = %s$ | quality: %s", rho_str, qual_str)
  } else {
    stat_val <- row_dt[["reg_status"]]
    stat_str <- if (!is.null(stat_val) && !is.na(stat_val)) as.character(stat_val) else "NA"
    subtitle_text <- sprintf("Spearman's $\\rho = %s$ | status: %s", rho_str, stat_str)
  }

  plot_data_list[[tw]] <- list(
    grid_dt = grid_dt,
    title_text = title_text,
    subtitle_text = subtitle_text,
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
# KROK 2: GENEROWANIE WYKRESÓW W JEDNYM PDF
##########################################################

pdf(pdf_file, width = 7, height = 7.5)

# --- CZĘŚĆ 1: 10 stron w TEJ SAMEJ SKALI (GLOBALNA) ---
for (tw in names(plot_data_list)) {
  item <- plot_data_list[[tw]]
  
  p <- ggtern(item$grid_dt, aes(x = x_A1, y = x_out, z = x_A2, color = density)) +
    geom_point(size = 1.5, stroke = 0) +
    labs(
      title = paste0(item$title_text, " (Global Scale)"),
      subtitle = item$subtitle_text,
      x = "p_A1",
      y = "p_out",
      z = "p_A2"
    ) +
    scale_L_continuous(breaks = item$breaks_LR / (1 - item$threshold), labels = item$labels_LR) +
    scale_T_continuous(breaks = (item$breaks_T - item$threshold) / (1 - item$threshold), labels = item$labels_T) +
    scale_R_continuous(breaks = item$breaks_LR / (1 - item$threshold), labels = item$labels_LR) +
    scale_color_viridis_c(
      option = "D", 
      name = density_legend_name, 
      limits = c(global_min_density, global_max_density),
      n.breaks = 5,
      guide = guide_colorbar(
        barwidth = unit(11, "cm"),
        barheight = unit(0.4, "cm"),
        title.position = "top",
        title.hjust = 0.5,
        label.theme = element_text(size = 7.5)
      )
    ) +
    theme_bw() +
    theme_showarrows() +
    theme_latex() +
    theme(
      tern.panel.expand = 0.6,
      tern.axis.arrow.sep = 0.4,
      legend.position = "bottom",
      legend.title = element_text(size = 9, face = "bold"),
      legend.margin = margin(t = 10),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5)
    )
    
  print(p)
}

# --- CZĘŚĆ 2: 10 stron w RÓŻNYCH SKALACH (LOKALNE) ---
for (tw in names(plot_data_list)) {
  item <- plot_data_list[[tw]]
  
  p <- ggtern(item$grid_dt, aes(x = x_A1, y = x_out, z = x_A2, color = density)) +
    geom_point(size = 1.5, stroke = 0) +
    labs(
      title = paste0(item$title_text, " (Local Scale)"),
      subtitle = item$subtitle_text,
      x = "p_A1",
      y = "p_out",
      z = "p_A2"
    ) +
    scale_L_continuous(breaks = item$breaks_LR / (1 - item$threshold), labels = item$labels_LR) +
    scale_T_continuous(breaks = (item$breaks_T - item$threshold) / (1 - item$threshold), labels = item$labels_T) +
    scale_R_continuous(breaks = item$breaks_LR / (1 - item$threshold), labels = item$labels_LR) +
    scale_color_viridis_c(
      option = "D", 
      name = density_legend_name,
      n.breaks = 5,
      guide = guide_colorbar(
        barwidth = unit(11, "cm"),
        barheight = unit(0.4, "cm"),
        title.position = "top",
        title.hjust = 0.5,
        label.theme = element_text(size = 7.5)
      )
    ) +
    theme_bw() +
    theme_showarrows() +
    theme_latex() +
    theme(
      tern.panel.expand = 0.6,
      tern.axis.arrow.sep = 0.4,
      legend.position = "bottom",
      legend.title = element_text(size = 9, face = "bold"),
      legend.margin = margin(t = 10),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5)
    )
    
  print(p)
}

dev.off()
cat(sprintf("Zapisano wielostronicowy plik PDF z poszerzoną legendą (%s): %s\n", scale_label, pdf_file))
