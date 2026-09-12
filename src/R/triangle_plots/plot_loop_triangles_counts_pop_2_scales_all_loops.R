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
  stop("Użycie: Rscript plot_loop_triangles_counts_pop_2_scales_all_loops.R <fit|reg> [log|log2|nolog]")
}

mode_type <- tolower(args[1]) # "fit" lub "reg"

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

# Zdefiniowane pętle do przetworzenia
target_loops <- c("L259", "L314", "L391", "L413")

log_str_info <- if (use_log) "LOGARYTMICZNA (log2)" else "LINIOWA"
cat(sprintf("Generowanie histogramów (2 skale) dla pętli [%s]: Tryb = %s | Skala = %s\n", 
            paste(target_loops, collapse = ", "), mode_type, log_str_info))

##########################################################
# KONFIGURACJA I ŚCIEŻKI
##########################################################

config_path <- file.path("config", "config.yml")
if (!file.exists(config_path)) {
  stop("Nie znaleziono pliku konfiguracji: ", config_path)
}
config <- yaml::yaml.load_file(config_path)
resultsdir <- config$paths$resultsdir

# Katalog ewaluacji z oknami czasowymi
eval_dir <- file.path(resultsdir, "MGLMreg", "eval_10h+_time_tissue_default")
if (!dir.exists(eval_dir)) {
  stop("Nie znaleziono katalogu ewaluacji: ", eval_dir)
}

# Pobranie wszystkich plików okien czasowych
eval_files <- list.files(eval_dir, pattern = "^eval_cor_.*\\.tsv\\.gz$", full.names = TRUE)

if (length(eval_files) == 0) {
  stop("Nie znaleziono żadnych plików ewaluacji w katalogu: ", eval_dir)
}

# Katalog docelowy PDF
out_dir <- file.path(resultsdir, "triangle_plots", "time_tissue", paste0("histograms_", mode_type))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

n_sim <- 100000
seed_base <- 12345
prefix <- paste0(mode_type, "_")

##########################################################
# GŁÓWNA PĘTLA PRZETWARZANIA
##########################################################

for (eval_file in eval_files) {
  # Wyciągnięcie nazwy okna czasowego z nazwy pliku
  tw_arg <- sub("^eval_cor_(.*)\\.tsv\\.gz$", "\\1", basename(eval_file))
  
  dt_eval <- fread(eval_file)
  dt_eval[, loop_id_str := as.character(loop_id)]
  dt_eval[, loop_id_clean := gsub("^L", "", loop_id_str)]
  
  for (loop_id_arg in target_loops) {
    dt_sub <- dt_eval[loop_id_clean == loop_id_arg | loop_id_str == loop_id_arg]
    
    if (nrow(dt_sub) == 0) {
      cat(sprintf("Brak danych dla loop_id=%s w oknie %s — pomijam.\n", loop_id_arg, tw_arg))
      next
    }
    
    cat(sprintf("Przetwarzanie: Loop = %s | TW = %s (%d populacji)\n", loop_id_arg, tw_arg, nrow(dt_sub)))
    setorder(dt_sub, population)
    
    pdf_filename <- if (use_log) {
      paste0("loop_", loop_id_arg, "_tw_", tw_arg, "_2scales_log2.pdf")
    } else {
      paste0("loop_", loop_id_arg, "_tw_", tw_arg, "_2scales.pdf")
    }
    pdf_file <- file.path(out_dir, pdf_filename)
    
    # --- FAZA 1: SYMULACJE, BUDOWANIE HISTOGRAMÓW I WYZNACZENIE SKALI GLOBALNEJ ---
    prepared_data <- list()
    all_counts <- c()
    
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
        cat(sprintf("  - Pominięto populację %s (brak parametrów dla x_out)\n", pop_id))
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
        cat(sprintf("  - Pominięto populację %s (brak parametrów drugiego kroku)\n", pop_id))
        next
      }
      
      # Symulacja z rozkładu GDM
      set.seed(seed_base + as.integer(gsub("[^0-9]", "", loop_id_arg)) + i)
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
      
      title_base <- sprintf("loop: %s | tw: %s | pop: %s (%s)", loop_id_arg, tw_arg, pop_id, t_name)
      
      if (mode_type == "fit") {
        qual_val <- row_dt[["fit_min_est_over_se"]]
        qual_str <- if (!is.null(qual_val) && !is.na(qual_val)) sprintf("%.3f", qual_val) else "NA"
        subtitle_text <- sprintf("Spearman's $\\rho = %s$ | quality: %s", rho_str, qual_str)
      } else {
        stat_val <- row_dt[["reg_status"]]
        stat_str <- if (!is.null(stat_val) && !is.na(stat_val)) as.character(stat_val) else "NA"
        subtitle_text <- sprintf("Spearman's $\\rho = %s$ | status: %s", rho_str, stat_str)
      }
      
      # Pomocnicza budowa wykresu w celu wyciągnięcia zliczeń (count) z geomu
      p_tmp <- ggtern(dts, aes(x = x_A1, y = x_out, z = x_A2)) +
        geom_tri_tern(aes(fill = after_stat(count)))
      
      gb <- ggplot_build(p_tmp)
      layer_counts <- gb$data[[1]]$count
      layer_counts <- layer_counts[!is.na(layer_counts)]
      
      all_counts <- c(all_counts, layer_counts)
      
      prepared_data[[pop_id]] <- list(
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
    
    if (length(prepared_data) == 0) {
      cat(sprintf("Brak poprawnych populacji do wygenerowania histogramów dla loop=%s, tw=%s\n", loop_id_arg, tw_arg))
      next
    }
    
    global_min <- min(all_counts, na.rm = TRUE)
    global_max <- max(all_counts, na.rm = TRUE)
    
    if (use_log) {
      pos_counts <- all_counts[all_counts > 0]
      global_min <- if (length(pos_counts) > 0) min(pos_counts, na.rm = TRUE) else 1
    }
    
    global_limits <- c(global_min, global_max)
    
    # --- FAZA 2: GENEROWANIE WYKRESÓW W PDF ---
    pdf(pdf_file, width = 6, height = 6.5)
    
    # CZĘŚĆ 1: SKALA GLOBALNA
    for (pop_id in names(prepared_data)) {
      item <- prepared_data[[pop_id]]
      
      scale_fill_global <- if (use_log) {
        scale_fill_viridis_c(option = "D", name = "Count (log2)", trans = "log2", limits = global_limits)
      } else {
        scale_fill_viridis_c(option = "D", name = "Count", limits = global_limits)
      }
      
      p_global <- ggtern(item$dts, aes(x = x_A1, y = x_out, z = x_A2)) +
        geom_tri_tern(aes(fill = after_stat(count))) +
        labs(
          title = paste0(item$title_base, " [GLOBAL SCALE]"),
          subtitle = item$subtitle_text,
          x = "p_A1",
          y = "p_out",
          z = "p_A2"
        ) +
        scale_L_continuous(breaks = item$breaks_LR / (1 - item$threshold), labels = item$labels_LR) +
        scale_T_continuous(breaks = (item$breaks_T - item$threshold) / (1 - item$threshold), labels = item$labels_T) +
        scale_R_continuous(breaks = item$breaks_LR / (1 - item$threshold), labels = item$labels_LR) +
        scale_fill_global +
        theme_bw() +
        theme_showarrows() +
        theme_latex() +
        theme(
          tern.panel.expand = 0.6,
          tern.axis.arrow.sep = 0.4,
          legend.position = "bottom",
          plot.title = element_text(hjust = 0.5, face = "bold", size = 11),
          plot.subtitle = element_text(hjust = 0.5, size = 10)
        )
      
      print(p_global)
    }
    
    # CZĘŚĆ 2: SKALA LOKALNA
    for (pop_id in names(prepared_data)) {
      item <- prepared_data[[pop_id]]
      
      scale_fill_local <- if (use_log) {
        scale_fill_viridis_c(option = "D", name = "Count (log2)", trans = "log2")
      } else {
        scale_fill_viridis_c(option = "D", name = "Count")
      }
      
      p_local <- ggtern(item$dts, aes(x = x_A1, y = x_out, z = x_A2)) +
        geom_tri_tern(aes(fill = after_stat(count))) +
        labs(
          title = paste0(item$title_base, " [LOCAL SCALE]"),
          subtitle = item$subtitle_text,
          x = "p_A1",
          y = "p_out",
          z = "p_A2"
        ) +
        scale_L_continuous(breaks = item$breaks_LR / (1 - item$threshold), labels = item$labels_LR) +
        scale_T_continuous(breaks = (item$breaks_T - item$threshold) / (1 - item$threshold), labels = item$labels_T) +
        scale_R_continuous(breaks = item$breaks_LR / (1 - item$threshold), labels = item$labels_LR) +
        scale_fill_local +
        theme_bw() +
        theme_showarrows() +
        theme_latex() +
        theme(
          tern.panel.expand = 0.6,
          tern.axis.arrow.sep = 0.4,
          legend.position = "bottom",
          plot.title = element_text(hjust = 0.5, face = "bold", size = 11),
          plot.subtitle = element_text(hjust = 0.5, size = 10)
        )
      
      print(p_local)
    }
    
    dev.off()
    cat(sprintf("Zapisano: %s\n", pdf_file))
  }
}

cat("Przetwarzanie zakończone pomyślnie dla wszystkich zdefiniowanych pętli i okien czasowych.\n")
