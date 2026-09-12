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

mode_type <- if (length(args) >= 1) tolower(args[1]) else "fit" # "fit" lub "reg"

if (!mode_type %in% c("fit", "reg")) {
  stop("Parametr 'mode_type' musi przyjmować wartość 'fit' lub 'reg'.")
}

# Obsługa argumentu 'log' / 'log2' (drugi argument, opcjonalny)
use_log <- FALSE
if (length(args) >= 2) {
  log_arg <- tolower(args[2])
  if (log_arg %in% c("log", "log2", "true", "1", "t", "yes")) {
    use_log <- TRUE
  }
}

log_str_info <- if (use_log) "LOGARYTMICZNA (log2)" else "LINIOWA"
cat(sprintf("Uruchamianie generowania gęstości dla wszystkich pętli | Tryb = %s | Skala = %s\n", mode_type, log_str_info))

##########################################################
# KONFIGURACJA I ŚCIEŻKI
##########################################################

config_path <- file.path("config", "config.yml")
if (!file.exists(config_path)) {
  stop("Nie znaleziono pliku konfiguracji: ", config_path)
}
config <- yaml::yaml.load_file(config_path)

resultsdir <- config$paths$resultsdir
datadir    <- config$paths$datadir

loops_file <- file.path(datadir, "long_and_short_range_loops_D_mel.tsv")
eval_dir   <- file.path(resultsdir, "MGLMreg", "eval_0h+_time_default")
out_dir    <- file.path(resultsdir, "triangle_plots", "time", paste0("density_", mode_type, "_batch"))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

time_windows <- c("00-02", "02-04", "04-06", "06-08", "08-10", 
                  "10-12", "12-14", "14-16", "16-18", "18-20")

prefix <- paste0(mode_type, "_")
grid_res <- 150
loops_per_pdf <- 10

##########################################################
# FUNKCJE POMOCNICZE
##########################################################

get_val_str <- function(dt, col_name, fmt = "%s") {
  if (!is.null(dt) && col_name %in% names(dt)) {
    val <- dt[[col_name]][1]
    if (!is.null(val) && !is.na(val) && length(val) > 0) {
      if (is.numeric(val)) {
        return(sprintf(fmt, val))
      }
      return(as.character(val))
    }
  }
  return("NA")
}

format_seconds <- function(seconds) {
  if (is.na(seconds) || is.infinite(seconds) || seconds < 0) return("--:--:--")
  h <- floor(seconds / 3600)
  m <- floor((seconds %% 3600) / 60)
  s <- floor(seconds %% 60)
  sprintf("%02d:%02d:%02d", h, m, s)
}

make_progress_bar <- function(percent, width = 20) {
  filled <- round((percent / 100) * width)
  empty <- width - filled
  paste0("[", paste(rep("=", filled), collapse = ""), paste(rep("-", empty), collapse = ""), "]")
}

dgdm3 <- function(x1, x2, x3, alpha, beta) {
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

  exp(log_pdf)
}

##########################################################
# 1. WCZYTANIE LISTY PĘTLI I PODZIAŁ NA PORCJE
##########################################################

if (!file.exists(loops_file)) {
  stop("Nie znaleziono pliku pętli: ", loops_file)
}

dt_loops <- fread(loops_file)
if (!"loop_id" %in% names(dt_loops)) {
  stop("Brak kolumny 'loop_id' w pliku pętli.")
}

all_loops <- unique(as.character(dt_loops$loop_id))
#all_loops <- all_loops[0:45]
all_loops <- all_loops[order(as.integer(sub("^L", "", all_loops)))]

num_loops <- length(all_loops)
cat(sprintf("Wczytano %d pętli z pliku.\n", num_loops))

loop_chunks <- split(all_loops, ceiling(seq_along(all_loops) / loops_per_pdf))
total_chunks <- length(loop_chunks)

##########################################################
# 2. WSTĘPNE WCZYTANIE DANYCH EVAL
##########################################################

cat("Wczytywanie danych ewaluacji do pamięci...\n")
eval_tables <- list()
for (tw in time_windows) {
  eval_file <- file.path(eval_dir, paste0("eval_cor_", tw, ".tsv.gz"))
  if (file.exists(eval_file)) {
    dt_tw <- fread(eval_file)
    dt_tw[, loop_id := as.character(loop_id)]
    eval_tables[[tw]] <- dt_tw
  } else {
    warning(sprintf("Pominięto okno %s (brak pliku %s)", tw, eval_file))
  }
}

##########################################################
# 3. GENEROWANIE PLIKÓW PDF DLA PORCJI PĘTLI
##########################################################

global_start_time <- Sys.time()
processed_loops_counter <- 0

for (chunk_idx in seq_along(loop_chunks)) {
  current_loops <- loop_chunks[[chunk_idx]]
  
  pdf_filename <- if (use_log) {
    sprintf("loops_part_%02d_density_%s_log2.pdf", chunk_idx, mode_type)
  } else {
    sprintf("loops_part_%02d_density_%s.pdf", chunk_idx, mode_type)
  }
  pdf_file <- file.path(out_dir, pdf_filename)
  
  cat(sprintf("\n========================================================================\n"))
  cat(sprintf(" ROZPOCZĘCIE PLIKU PDF %d/%d (%d pętli) -> %s\n", 
              chunk_idx, total_chunks, length(current_loops), basename(pdf_file)))
  cat(sprintf("========================================================================\n"))
  
  pdf(pdf_file, width = 7, height = 7.5)
  
  for (l_idx in seq_along(current_loops)) {
    loop_id_arg <- current_loops[l_idx]
    processed_loops_counter <- processed_loops_counter + 1
    
    elapsed_secs <- as.numeric(difftime(Sys.time(), global_start_time, units = "secs"))
    pct <- (processed_loops_counter / num_loops) * 100
    avg_secs_per_loop <- elapsed_secs / processed_loops_counter
    remaining_loops <- num_loops - processed_loops_counter
    eta_secs <- remaining_loops * avg_secs_per_loop
    
    pbar <- make_progress_bar(pct)
    
    cat(sprintf("%s %5.1f%% | Pętla %4d/%d [%s] | Plik %d/%d (%2d/%2d) | Upłynęło: %s | ETA: %s\n",
                pbar,
                pct,
                processed_loops_counter,
                num_loops,
                loop_id_arg,
                chunk_idx,
                total_chunks,
                l_idx,
                length(current_loops),
                format_seconds(elapsed_secs),
                format_seconds(eta_secs)))
    
    plot_data_list <- list()
    global_min_density <- Inf
    global_max_density <- -Inf
    
    for (tw in time_windows) {
      if (!tw %in% names(eval_tables)) next
      dt_eval <- eval_tables[[tw]]
      row_dt <- dt_eval[loop_id == loop_id_arg]
      if (nrow(row_dt) == 0) next
      row_dt <- row_dt[1]
      
      # Pobranie parametrów opisu niezależnie od powodzenia estymacji
      rho_str <- get_val_str(row_dt, paste0(prefix, "spearman_rho"), "%.3f")
      stat_str <- get_val_str(row_dt, "reg_status")
      if (stat_str == "NA") stat_str <- get_val_str(row_dt, "status")
      qual_str <- get_val_str(row_dt, "fit_min_est_over_se", "%.3f")
      
      a_out <- row_dt[[paste0(prefix, "alpha_x_out_est")]]
      b_out <- row_dt[[paste0(prefix, "beta_x_out_est")]]
      a_A1  <- row_dt[[paste0(prefix, "alpha_x_A1_est")]]
      b_A1  <- row_dt[[paste0(prefix, "beta_x_A1_est")]]
      a_A2  <- row_dt[[paste0(prefix, "alpha_x_A2_est")]]
      b_A2  <- row_dt[[paste0(prefix, "beta_x_A2_est")]]
      
      # Weryfikacja poprawności estymacji parametru
      if (is.null(a_out) || is.null(b_out) || is.na(a_out) || is.na(b_out) || a_out <= 0 || b_out <= 0) next
      
      if (!is.null(a_A1) && !is.null(b_A1) && !is.na(a_A1) && !is.na(b_A1) && a_A1 > 0 && b_A1 > 0) {
        cat2_name <- "x_A1"; cat3_name <- "x_A2"; a2 <- a_A1; b2 <- b_A1
      } else if (!is.null(a_A2) && !is.null(b_A2) && !is.na(a_A2) && !is.na(b_A2) && a_A2 > 0 && b_A2 > 0) {
        cat2_name <- "x_A2"; cat3_name <- "x_A1"; a2 <- a_A2; b2 <- b_A2
      } else {
        next
      }
      
      alpha_vec <- c(a_out, a2)
      beta_vec  <- c(b_out, b2)
      mean_out  <- a_out / (a_out + b_out)
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
      
      grid_dt[, density := dgdm3(val_1, val_2, val_3, alpha_vec, beta_vec)]
      
      min_d <- min(grid_dt$density, na.rm = TRUE)
      max_d <- max(grid_dt$density, na.rm = TRUE)
      
      if (use_log) {
        pos_densities <- grid_dt$density[grid_dt$density > 0]
        if (length(pos_densities) > 0) min_d <- min(pos_densities, na.rm = TRUE) else min_d <- 1e-12
      }
      
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
      
      title_text <- sprintf("loop_id: %s | time_window: %s", loop_id_arg, tw)
      
      if (mode_type == "fit") {
        subtitle_text <- sprintf("Spearman's $\\rho = %s$ | quality: %s", rho_str, qual_str)
      } else {
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
    
    if (length(plot_data_list) == 0) next
    
    color_title <- if (use_log) "Density (log2)" else "Density"
    
    # --- Rysowanie: Skala Globalna ---
    for (tw in names(plot_data_list)) {
      item <- plot_data_list[[tw]]
      
      color_scale_global <- if (use_log) {
        scale_color_viridis_c(
          option = "D", 
          name = color_title,
          trans = "log2",
          limits = c(global_min_density, global_max_density),
          n.breaks = 5,
          guide = guide_colorbar(
            barwidth = unit(11, "cm"),
            barheight = unit(0.4, "cm"),
            title.position = "top",
            title.hjust = 0.5,
            label.theme = element_text(size = 7.5)
          )
        )
      } else {
        scale_color_viridis_c(
          option = "D", 
          name = color_title,
          limits = c(global_min_density, global_max_density),
          n.breaks = 5,
          guide = guide_colorbar(
            barwidth = unit(11, "cm"),
            barheight = unit(0.4, "cm"),
            title.position = "top",
            title.hjust = 0.5,
            label.theme = element_text(size = 7.5)
          )
        )
      }
      
      p <- ggtern(item$grid_dt, aes(x = x_A1, y = x_out, z = x_A2, color = density)) +
        geom_point(size = 1.5, stroke = 0) +
        labs(
          title = paste0(item$title_text, " (Global Scale)"),
          subtitle = item$subtitle_text,
          x = "p_A1", y = "p_out", z = "p_A2"
        ) +
        scale_L_continuous(breaks = item$breaks_LR / (1 - item$threshold), labels = item$labels_LR) +
        scale_T_continuous(breaks = (item$breaks_T - item$threshold) / (1 - item$threshold), labels = item$labels_T) +
        scale_R_continuous(breaks = item$breaks_LR / (1 - item$threshold), labels = item$labels_LR) +
        color_scale_global +
        theme_bw() + theme_showarrows() + theme_latex() +
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
    
    # --- Rysowanie: Skala Lokalna ---
    for (tw in names(plot_data_list)) {
      item <- plot_data_list[[tw]]
      
      color_scale_local <- if (use_log) {
        scale_color_viridis_c(
          option = "D", 
          name = color_title,
          trans = "log2",
          n.breaks = 5,
          guide = guide_colorbar(
            barwidth = unit(11, "cm"),
            barheight = unit(0.4, "cm"),
            title.position = "top",
            title.hjust = 0.5,
            label.theme = element_text(size = 7.5)
          )
        )
      } else {
        scale_color_viridis_c(
          option = "D", 
          name = color_title,
          n.breaks = 5,
          guide = guide_colorbar(
            barwidth = unit(11, "cm"),
            barheight = unit(0.4, "cm"),
            title.position = "top",
            title.hjust = 0.5,
            label.theme = element_text(size = 7.5)
          )
        )
      }
      
      p <- ggtern(item$grid_dt, aes(x = x_A1, y = x_out, z = x_A2, color = density)) +
        geom_point(size = 1.5, stroke = 0) +
        labs(
          title = paste0(item$title_text, " (Local Scale)"),
          subtitle = item$subtitle_text,
          x = "p_A1", y = "p_out", z = "p_A2"
        ) +
        scale_L_continuous(breaks = item$breaks_LR / (1 - item$threshold), labels = item$labels_LR) +
        scale_T_continuous(breaks = (item$breaks_T - item$threshold) / (1 - item$threshold), labels = item$labels_T) +
        scale_R_continuous(breaks = item$breaks_LR / (1 - item$threshold), labels = item$labels_LR) +
        color_scale_local +
        theme_bw() + theme_showarrows() + theme_latex() +
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
  }
  
  dev.off()
  cat(sprintf("-> Zakończono i zapisano plik PDF: %s\n", pdf_file))
}

total_elapsed <- as.numeric(difftime(Sys.time(), global_start_time, units = "secs"))
cat(sprintf("\n========================================================================\n"))
cat(sprintf(" ZAKOŃCZONO GENEROWANIE WSZYSTKICH WYKRESÓW!\n"))
cat(sprintf(" Przetworzono pętli: %d | Całkowity czas: %s\n", num_loops, format_seconds(total_elapsed)))
cat(sprintf("========================================================================\n"))
