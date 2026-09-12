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
  stop("Użycie: Rscript plot_loop_triangles_counts.R <loop_id> <fit|reg> [log|nolog]")
}

loop_id_arg <- as.character(args[1])
mode_type   <- tolower(args[2]) # "fit" lub "reg"

if (!mode_type %in% c("fit", "reg")) {
  stop("Parametr 'mode_type' musi przyjmować wartość 'fit' lub 'reg'.")
}

# Obsługa argumentu 'log' (trzeci argument, opcjonalny)
use_log <- FALSE
if (length(args) >= 3) {
  log_arg <- tolower(args[3])
  if (log_arg %in% c("log", "log2", "true", "1", "t", "yes")) {
    use_log <- TRUE
  }
}

log_str_info <- if (use_log) "LOGARYTMICZNA (log2)" else "LINIOWA"
cat(sprintf("Uruchamianie generowania trójkątów z plików eval: Loop = %s | Tryb = %s | Skala = %s\n", 
            loop_id_arg, mode_type, log_str_info))

##########################################################
# KONFIGURACJA I ŚCIEŻKI
##########################################################

config_path <- file.path("config", "config.yml")
if (!file.exists(config_path)) {
  stop("Nie znaleziono pliku konfiguracji: ", config_path)
}
config <- yaml::yaml.load_file(config_path)
resultsdir <- config$paths$resultsdir

# Katalog i plik docelowy PDF (dodanie sufiksu _log2 w przypadku skali logarytmicznej)
out_dir <- file.path(resultsdir, "triangle_plots", "time", paste0("histograms_", mode_type))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

pdf_filename <- if (use_log) paste0("loop_", loop_id_arg, "_log2.pdf") else paste0("loop_", loop_id_arg, ".pdf")
pdf_file <- file.path(out_dir, pdf_filename)

# Ścieżka do katalogu z wynikami ewaluacji okien czasowych
eval_dir <- file.path(resultsdir, "MGLMreg", "eval_0h+_time_default")

time_windows <- c("00-02", "02-04", "04-06", "06-08", "08-10", 
                  "10-12", "12-14", "14-16", "16-18", "18-20")

n_sim <- 100000
seed_base <- 12345
prefix <- paste0(mode_type, "_")

##########################################################
# GENEROWANIE WYKRESÓW W PDF
##########################################################

pdf(pdf_file, width = 6, height = 6.5)

for (tw in time_windows) {
  eval_file <- file.path(eval_dir, paste0("eval_cor_", tw, ".tsv.gz"))
  
  if (!file.exists(eval_file)) {
    cat(sprintf("Pominięto okno %s (brak pliku eval: %s)\n", tw, eval_file))
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
  
  # Ekstrakcja parametrów rozkładu Beta z pliku ewaluacji
  a_out <- row_dt[[paste0(prefix, "alpha_x_out_est")]]
  b_out <- row_dt[[paste0(prefix, "beta_x_out_est")]]
  
  a_A1  <- row_dt[[paste0(prefix, "alpha_x_A1_est")]]
  b_A1  <- row_dt[[paste0(prefix, "beta_x_A1_est")]]
  
  a_A2  <- row_dt[[paste0(prefix, "alpha_x_A2_est")]]
  b_A2  <- row_dt[[paste0(prefix, "beta_x_A2_est")]]
  
  if (is.na(a_out) || is.na(b_out) || a_out <= 0 || b_out <= 0) {
    cat(sprintf("Pominięto okno %s (brak poprawnych estymat parametrów dla x_out)\n", tw))
    next
  }
  
  # Określenie kolejności w kaskadzie GDM
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
    cat(sprintf("Pominięto okno %s (brak parametrów dla drugiego kroku kaskady)\n", tw))
    next
  }
  
  # Symulacja punktów z rozkładu GDM
  set.seed(seed_base + as.integer(gsub("[^0-9]", "", loop_id_arg)))
  p1 <- rbeta(n_sim, a_out, b_out)
  p2 <- rbeta(n_sim, a2, b2) * (1 - p1)
  p3 <- pmax(0, 1 - p1 - p2)
  
  dt_sim <- data.table(x_out = p1)
  dt_sim[, (cat2_name) := p2]
  dt_sim[, (cat3_name) := p3]
  
  # Przeskalowanie trójkąta
  mean_out <- mean(dt_sim$x_out)
  threshold <- 2 * mean_out - 1
  
  dts <- dt_sim[x_out >= threshold, ]
  dts[, x_out := (x_out - threshold) / (1 - threshold)]
  dts[, x_A1  := x_A1 / (1 - threshold)]
  dts[, x_A2  := x_A2 / (1 - threshold)]
  
  # Podziałki i etykiety osi
  breaks_LR <- pretty(c(0, 1 - threshold))
  breaks_LR <- breaks_LR[breaks_LR > 0 & breaks_LR <= 1 - threshold]
  labels_LR <- sprintf("%g", breaks_LR)
  labels_LR <- paste0("$", sub("-0", "-", sub("e(.*)", " \\\\cdot 10^{\\1}", labels_LR)), "$")

  breaks_T  <- pretty(c(threshold, 1))
  breaks_T  <- breaks_T[breaks_T >= threshold]
  labels_T  <- sub(" - 0$", "", sprintf("1 - %g", 1 - breaks_T))
  labels_T  <- paste0("$", sub("-0", "-", sub("e(.*)", " \\\\cdot 10^{\\1}", labels_T)), "$")
  
  # Odczyt informacji do nagłówków z pliku eval
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
  
  # Konfiguracja skali wypełnienia (liniowa vs logarytmiczna log2)
  fill_scale <- if (use_log) {
    scale_fill_viridis_c(option = "D", name = "Count (log2)", trans = "log2")
  } else {
    scale_fill_viridis_c(option = "D", name = "Count")
  }

  # Tworzenie wykresu
  p <- ggtern(dts, aes(x = x_A1, y = x_out, z = x_A2)) +
    geom_tri_tern(aes(fill = after_stat(count))) +
    labs(
      title = title_text,
      subtitle = subtitle_text,
      x = "p_A1",
      y = "p_out",
      z = "p_A2"
    ) +
    scale_L_continuous(breaks = breaks_LR / (1 - threshold), labels = labels_LR) +
    scale_T_continuous(breaks = (breaks_T - threshold) / (1 - threshold), labels = labels_T) +
    scale_R_continuous(breaks = breaks_LR / (1 - threshold), labels = labels_LR) +
    fill_scale +
    theme_bw() +
    theme_showarrows() +
    theme_latex() +
    theme(
      tern.panel.expand = 0.6,
      tern.axis.arrow.sep = 0.4,
      legend.position = "bottom",
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5)
    )
    
  print(p)
}

dev.off()
cat(sprintf("Zapisano wielostronicowy plik PDF: %s\n", pdf_file))
