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

if (length(args) < 1) {
  stop("Użycie: Rscript plot_loop_triangles_counts_2_scales_all_loops.R <fit|reg> [log|log2|nolog] [loops_per_pdf]")
}

mode_type <- tolower(args[1]) # "fit" lub "reg"

if (!mode_type %in% c("fit", "reg")) {
  stop("Parametr 'mode_type' musi przyjmować wartość 'fit' lub 'reg'.")
}

# Obsługa argumentu skali (trzeci argument)
use_log <- FALSE
if (length(args) >= 2) {
  log_arg <- tolower(args[2])
  if (log_arg %in% c("log", "log2", "true", "1", "t", "yes")) {
    use_log <- TRUE
  }
}

# Liczba pętli na jeden plik PDF (domyślnie 80)
loops_per_pdf <- 80
if (length(args) >= 3) {
  parsed_lpp <- as.integer(args[3])
  if (!is.na(parsed_lpp) && parsed_lpp > 0) {
    loops_per_pdf <- parsed_lpp
  }
}

log_str_info <- if (use_log) "LOGARYTMICZNA (log2)" else "LINIOWA"
cat(sprintf("Uruchamianie przetwarzania wszystkich pętli (2 skale)\nTryb: %s | Skala: %s | Pętli na PDF: %d\n", 
            mode_type, log_str_info, loops_per_pdf))

##########################################################
# KONFIGURACJA I ŚCIEŻKI
##########################################################

config_path <- file.path("config", "config.yml")
if (!file.exists(config_path)) {
  stop("Nie znaleziono pliku konfiguracji: ", config_path)
}
config <- yaml::yaml.load_file(config_path)
resultsdir <- config$paths$resultsdir

out_dir <- file.path(resultsdir, "triangle_plots", "time", paste0("histograms_", mode_type, "_all"))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

eval_dir <- file.path(resultsdir, "MGLMreg", "eval_0h+_time_default")

time_windows <- c("00-02", "02-04", "04-06", "06-08", "08-10", 
                  "10-12", "12-14", "14-16", "16-18", "18-20")

n_sim <- 100000
seed_base <- 12345
prefix <- paste0(mode_type, "_")

##########################################################
# POBLIERANIE LISTY WSZYSTKICH UNIKALNYCH LOOP_ID
##########################################################

cat("Wyszukiwanie unikalnych loop_id...\n")
all_loops <- character()

for (tw in time_windows) {
  eval_file <- file.path(eval_dir, paste0("eval_cor_", tw, ".tsv.gz"))
  if (file.exists(eval_file)) {
    dt_tmp <- fread(eval_file, select = "loop_id")
    all_loops <- unique(c(all_loops, as.character(dt_tmp$loop_id)))
  }
}

all_loops <- sort(all_loops)
total_loops <- length(all_loops)

if (total_loops == 0) {
  stop("Nie znaleziono żadnych pętli w plikach ewaluacyjnych.")
}

cat(sprintf("Znaleziono łącznie %d pętli.\n", total_loops))

# Podział pętli na pakiety (batches)
loop_batches <- split(all_loops, ceiling(seq_along(all_loops) / loops_per_pdf))
total_batches <- length(loop_batches)

##########################################################
# PĘTLA GŁÓWNA - GENEROWANIE PLIKÓW PDF (BATCHE)
##########################################################

for (b_idx in seq_along(loop_batches)) {
  current_batch_loops <- loop_batches[[b_idx]]
  
  log_suffix <- if (use_log) "_log2" else ""
  pdf_filename <- sprintf("loops_batch_%02d_of_%02d%s.pdf", b_idx, total_batches, log_suffix)
  pdf_file <- file.path(out_dir, pdf_filename)
  
  cat(sprintf("\n--- Przetwarzanie pakietu %d/%d (%d pętli) -> %s ---\n", 
              b_idx, total_batches, length(current_batch_loops), pdf_filename))
  
  pdf(pdf_file, width = 6, height = 6.5)
  
  for (l_idx in seq_along(current_batch_loops)) {
    loop_id_arg <- current_batch_loops[l_idx]
    
    # --------------------------------------------------------
    # PRZYGOTOWANIE DANYCH DLA DANEJ PĘTLI
    # --------------------------------------------------------
    prepared_data <- list()
    all_counts <- c()
    
    for (tw in time_windows) {
      eval_file <- file.path(eval_dir, paste0("eval_cor_", tw, ".tsv.gz"))
      if (!file.exists(eval_file)) next
      
      dt_eval <- fread(eval_file)
      dt_eval[, loop_id := as.character(loop_id)]
      row_dt <- dt_eval[loop_id == loop_id_arg]
      if (nrow(row_dt) == 0) next
      row_dt <- row_dt[1]
      
      a_out <- row_dt[[paste0(prefix, "alpha_x_out_est")]]
      b_out <- row_dt[[paste0(prefix, "beta_x_out_est")]]
      a_A1  <- row_dt[[paste0(prefix, "alpha_x_A1_est")]]
      b_A1  <- row_dt[[paste0(prefix, "beta_x_A1_est")]]
      a_A2  <- row_dt[[paste0(prefix, "alpha_x_A2_est")]]
      b_A2  <- row_dt[[paste0(prefix, "beta_x_A2_est")]]
      
      if (is.na(a_out) || is.na(b_out) || a_out <= 0 || b_out <= 0) next
      
      if (!is.na(a_A1) && !is.na(b_A1) && a_A1 > 0 && b_A1 > 0) {
        cat2_name <- "x_A1"; cat3_name <- "x_A2"
        a2 <- a_A1; b2 <- b_A1
      } else if (!is.na(a_A2) && !is.na(b_A2) && a_A2 > 0 && b_A2 > 0) {
        cat2_name <- "x_A2"; cat3_name <- "x_A1"
        a2 <- a_A2; b2 <- b_A2
      } else {
        next
      }
      
      set.seed(seed_base + as.integer(gsub("[^0-9]", "", loop_id_arg)))
      p1 <- rbeta(n_sim, a_out, b_out)
      p2 <- rbeta(n_sim, a2, b2) * (1 - p1)
      p3 <- pmax(0, 1 - p1 - p2)
      
      dt_sim <- data.table(x_out = p1)
      dt_sim[, (cat2_name) := p2]
      dt_sim[, (cat3_name) := p3]
      
      mean_out <- mean(dt_sim$x_out)
      threshold <- 2 * mean_out - 1
      
      dts <- dt_sim[x_out >= threshold, ]
      dts[, x_out := (x_out - threshold) / (1 - threshold)]
      dts[, x_A1  := x_A1 / (1 - threshold)]
      dts[, x_A2  := x_A2 / (1 - threshold)]
      
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
      
      title_base <- sprintf("loop_id: %s | time_window: %s", loop_id_arg, tw)
      
      if (mode_type == "fit") {
        qual_val <- row_dt[["fit_min_est_over_se"]]
        qual_str <- if (!is.null(qual_val) && !is.na(qual_val)) sprintf("%.3f", qual_val) else "NA"
        subtitle_text <- sprintf("Spearman's $\\rho = %s$ | quality: %s", rho_str, qual_str)
      } else {
        stat_val <- row_dt[["reg_status"]]
        stat_str <- if (!is.null(stat_val) && !is.na(stat_val)) as.character(stat_val) else "NA"
        subtitle_text <- sprintf("Spearman's $\\rho = %s$ | status: %s", rho_str, stat_str)
      }
      
      # Wyznaczenie liczności binów histogramu dla skali globalnej danej pętli
      p_temp <- ggtern(dts, aes(x = x_A1, y = x_out, z = x_A2)) +
        geom_tri_tern(aes(fill = after_stat(count)))
      
      gb <- ggplot_build(p_temp)
      counts <- gb$data[[1]]$count
      all_counts <- c(all_counts, counts[!is.na(counts)])
      
      prepared_data[[tw]] <- list(
        dts = dts,
        breaks_LR = breaks_LR,
        labels_LR = labels_LR,
        breaks_T = breaks_T,
        labels_T = labels_T,
        threshold = threshold,
        title_base = title_base,
        subtitle_text = subtitle_text
      )
    }
    
    if (length(prepared_data) == 0) next
    
    # Wyznaczenie granic skali globalnej dla bieżącej pętli
    global_min <- min(all_counts, na.rm = TRUE)
    global_max <- max(all_counts, na.rm = TRUE)
    if (use_log) global_min <- max(1, global_min)
    global_limits <- c(global_min, global_max)
    
    # --------------------------------------------------------
    # RENDEROWANIE WYKRESÓW DLA PĘTLI
    # --------------------------------------------------------
    
    # 1. 10 okien w skali globalnej dla tej pętli
    for (tw in names(prepared_data)) {
      item <- prepared_data[[tw]]
      
      fill_scale_global <- if (use_log) {
        scale_fill_viridis_c(option = "D", name = "Count (log2)", trans = "log2", limits = global_limits)
      } else {
        scale_fill_viridis_c(option = "D", name = "Count", limits = global_limits)
      }
      
      p_global <- ggtern(item$dts, aes(x = x_A1, y = x_out, z = x_A2)) +
        geom_tri_tern(aes(fill = after_stat(count))) +
        labs(
          title = paste0(item$title_base, " [GLOBAL SCALE]"),
          subtitle = item$subtitle_text,
          x = "p_A1", y = "p_out", z = "p_A2"
        ) +
        scale_L_continuous(breaks = item$breaks_LR / (1 - item$threshold), labels = item$labels_LR) +
        scale_T_continuous(breaks = (item$breaks_T - item$threshold) / (1 - item$threshold), labels = item$labels_T) +
        scale_R_continuous(breaks = item$breaks_LR / (1 - item$threshold), labels = item$labels_LR) +
        fill_scale_global +
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
        
      print(p_global)
    }
    
    # 2. 10 okien w skalach lokalnych dla tej pętli
    for (tw in names(prepared_data)) {
      item <- prepared_data[[tw]]
      
      fill_scale_local <- if (use_log) {
        scale_fill_viridis_c(option = "D", name = "Count (log2)", trans = "log2")
      } else {
        scale_fill_viridis_c(option = "D", name = "Count")
      }
      
      p_local <- ggtern(item$dts, aes(x = x_A1, y = x_out, z = x_A2)) +
        geom_tri_tern(aes(fill = after_stat(count))) +
        labs(
          title = paste0(item$title_base, " [LOCAL SCALE]"),
          subtitle = item$subtitle_text,
          x = "p_A1", y = "p_out", z = "p_A2"
        ) +
        scale_L_continuous(breaks = item$breaks_LR / (1 - item$threshold), labels = item$labels_LR) +
        scale_T_continuous(breaks = (item$breaks_T - item$threshold) / (1 - item$threshold), labels = item$labels_T) +
        scale_R_continuous(breaks = item$breaks_LR / (1 - item$threshold), labels = item$labels_LR) +
        fill_scale_local +
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
        
      print(p_local)
    }
    
    if (l_idx %% 10 == 0 || l_idx == length(current_batch_loops)) {
      cat(sprintf("   Wygenerowano: %d/%d pętli w bieżącym pakiecie\n", l_idx, length(current_batch_loops)))
    }
  }
  
  dev.off()
  cat(sprintf("Zapisano plik: %s\n", pdf_file))
}

cat("\nPrzetwarzanie zakończone pomyślnie!\n")
