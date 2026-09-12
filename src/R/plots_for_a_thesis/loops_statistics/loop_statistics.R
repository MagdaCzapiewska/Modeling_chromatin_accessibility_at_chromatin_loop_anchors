#!/usr/bin/env Rscript

library(data.table)
library(ggplot2)
library(yaml)
library(scales)
library(patchwork)

# ==============================================================================
# 1. KONFIGURACJA I ŚCIEŻKI
# ==============================================================================
config_path <- file.path("config", "config.yml")
if (!file.exists(config_path)) {
  stop("Nie znaleziono pliku konfiguracji: ", config_path)
}
config <- yaml::yaml.load_file(config_path)

data_dir    <- config$paths$datadir
results_dir <- config$paths$resultsdir

output_dir <- file.path(results_dir, "plots_for_a_thesis", "loop_statistics")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

loops_file <- file.path(data_dir, "long_and_short_range_loops_D_mel.tsv")

output_violin               <- file.path(output_dir, "violin_loop_length_log10.pdf")
output_histogram_loops       <- file.path(output_dir, "histogram_loop_length.pdf")
output_histogram_anchors_full<- file.path(output_dir, "histogram_anchor_length.pdf")
output_histogram_anchors_20k <- file.path(output_dir, "histogram_anchor_length_0_20k.pdf")

# ==============================================================================
# 2. WCZYTYWANIE I PRZYGOTOWANIE DANYCH
# ==============================================================================
if (!file.exists(loops_file)) {
  stop("Nie znaleziono pliku z pętlami: ", loops_file)
}

cat("Wczytywanie danych z pliku:", loops_file, "\n")
dt <- fread(loops_file)

# Obliczanie długości kotwic (Anchor 1 i Anchor 2)
dt[, anchor1_len := abs(x2 - x1)]
dt[, anchor2_len := abs(y2 - y1)]

# Obliczanie długości pętli (odległość między środkami kotwic)
dt[, loop_length := abs((y1 + y2) / 2 - (x1 + x2) / 2)]
dt[, log10_loop_length := log10(loop_length)]

# Przygotowanie połączonej tabeli długości kotwic (wszystkie kotwice razem)
dt_anchors <- rbind(
  dt[, .(anchor_len = anchor1_len, anchor_type = "Anchor 1")],
  dt[, .(anchor_len = anchor2_len, anchor_type = "Anchor 2")]
)

# ==============================================================================
# 3. OBLICZANIE I WYPISYWANIE STATYSTYK NA KONSOLĘ
# ==============================================================================
mean_loop_len   <- mean(dt$loop_length, na.rm = TRUE)
median_loop_len <- median(dt$loop_length, na.rm = TRUE)

mean_anchor_len   <- mean(dt_anchors$anchor_len, na.rm = TRUE)
median_anchor_len <- median(dt_anchors$anchor_len, na.rm = TRUE)

cat("\n==============================================================================\n")
cat("                       STATYSTYKI DŁUGOŚCI PĘTLI I KOTWIC                    \n")
cat("==============================================================================\n")
cat(sprintf("Długość pętli (Loop Length):\n"))
cat(sprintf("  - Średnia (Mean):   %s bp\n", format(round(mean_loop_len, 2), big.mark = ",")))
cat(sprintf("  - Mediana (Median): %s bp\n", format(round(median_loop_len, 2), big.mark = ",")))
cat("\n")
cat(sprintf("Długość kotwicy (Anchor Length):\n"))
cat(sprintf("  - Średnia (Mean):   %s bp\n", format(round(mean_anchor_len, 2), big.mark = ",")))
cat(sprintf("  - Mediana (Median): %s bp\n", format(round(median_anchor_len, 2), big.mark = ",")))
cat("==============================================================================\n\n")

# ==============================================================================
# MOTYW GRAFICZNY (PUB-READY THEME)
# ==============================================================================
custom_theme <- theme_minimal(base_size = 12) +
  theme(
    plot.title       = element_text(size = 13, face = "bold", hjust = 0.5, color = "black", margin = margin(b = 6)),
    plot.subtitle    = element_text(size = 12.5, hjust = 0.5, color = "black", margin = margin(b = 8)),
    axis.title       = element_text(size = 11, face = "bold", color = "black"),
    axis.text        = element_text(size = 10, color = "black"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "#f1f2f6"),
    panel.border     = element_rect(color = "#dcdde1", fill = NA, linewidth = 0.5)
  )

# ==============================================================================
# WYKRES 1: WYKRES WIOLINOWY log10(DŁUGOŚĆ PĘTLI)
# ==============================================================================
cat("Generowanie wykresu wiolinowego...\n")

p_violin <- ggplot(dt, aes(x = "", y = log10_loop_length)) +
  geom_violin(fill = "#4682B4", color = "black", alpha = 0.7, linewidth = 0.5) +
  geom_boxplot(width = 0.12, fill = "white", color = "black", outlier.size = 0.8, outlier.alpha = 0.5) +
  scale_y_continuous(labels = label_comma()) +
  labs(
    title = expression(bold("Distribution of chromatin loop lengths (log"[10]*")")),
    subtitle = paste0("Total number of loops: ", format(nrow(dt), big.mark = ",")),
    x = NULL,
    y = expression(bold("log"[10]*" (loop length in bp)"))
  ) +
  custom_theme +
  theme(
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank()
  )

pdf(output_violin, width = 5, height = 6)
print(p_violin)
dev.off()

# ==============================================================================
# WYKRES 2: HISTOGRAM DŁUGOŚCI PĘTLI (SKALA NORMALNA I LOG10)
# ==============================================================================
cat("Generowanie histogramu długości pętli...\n")

p_hist_linear <- ggplot(dt, aes(x = loop_length)) +
  geom_histogram(bins = 40, fill = "#3498db", color = "white", linewidth = 0.2) +
  scale_x_continuous(labels = label_comma()) +
  scale_y_continuous(labels = label_comma()) +
  labs(
    title = "Linear Scale",
    x = "Loop length (bp)",
    y = "Number of loops in a bin"
  ) +
  custom_theme

p_hist_log10 <- ggplot(dt, aes(x = loop_length)) +
  geom_histogram(bins = 40, fill = "#2ecc71", color = "white", linewidth = 0.2) +
  scale_x_log10(labels = label_comma()) +
  scale_y_continuous(labels = label_comma()) +
  labs(
    title = expression("Log"[10]*" Scale"),
    x = "Loop length (bp, log scale)",
    y = "Number of loops in a bin"
  ) +
  custom_theme

combined_hist_loops <- (p_hist_linear | p_hist_log10) +
  plot_annotation(
    title = "Chromatin Loop Length Distribution",
    subtitle = paste0("Total number of loops: ", format(nrow(dt), big.mark = ",")),
    theme = theme(
      plot.title    = element_text(size = 14, face = "bold", hjust = 0.5, color = "black"),
      plot.subtitle = element_text(size = 12.5, hjust = 0.5, color = "black")
    )
  )

pdf(output_histogram_loops, width = 10, height = 4.5)
print(combined_hist_loops)
dev.off()

# ==============================================================================
# WYKRES 3A: HISTOGRAM DŁUGOŚCI KOTWIC (PEŁNY ZAKRES)
# ==============================================================================
cat("Generowanie histogramu długości kotwic (pełny zakres)...\n")

p_hist_anchors_full <- ggplot(dt_anchors, aes(x = anchor_len)) +
  geom_histogram(bins = 50, fill = "#e74c3c", color = "white", linewidth = 0.2) +
  scale_x_continuous(labels = label_comma()) +
  scale_y_continuous(labels = label_comma()) +
  geom_vline(xintercept = median_anchor_len, linetype = "dashed", color = "black", linewidth = 0.8) +
  annotate(
    "text", 
    x = median_anchor_len, 
    y = Inf, 
    label = paste0(" Overall median: ", format(round(median_anchor_len, 1), big.mark = ","), " bp"), 
    vjust = 2, 
    hjust = -0.05, 
    fontface = "bold", 
    color = "black"
  ) +
  labs(
    title = "Anchor length distribution",
    subtitle = paste0("Total number of anchors: ", format(nrow(dt_anchors), big.mark = ",")),
    x = "Anchor length (bp)",
    y = "Number of anchors in a bin"
  ) +
  custom_theme

pdf(output_histogram_anchors_full, width = 7, height = 5)
print(p_hist_anchors_full)
dev.off()

# ==============================================================================
# WYKRES 3B: HISTOGRAM DŁUGOŚCI KOTWIC (ZAKRES 0 – 20,000 bp)
# ==============================================================================
cat("Generowanie histogramu długości kotwic (zakres 0 - 20,000 bp)...\n")

dt_anchors_20k <- dt_anchors[anchor_len <= 20000]

p_hist_anchors_20k <- ggplot(dt_anchors_20k, aes(x = anchor_len)) +
  geom_histogram(bins = 50, fill = "#e74c3c", color = "white", linewidth = 0.2) +
  scale_x_continuous(limits = c(0, 20000), labels = label_comma()) +
  scale_y_continuous(labels = label_comma()) +
  geom_vline(xintercept = median_anchor_len, linetype = "dashed", color = "black", linewidth = 0.8) +
  annotate(
    "text", 
    x = median_anchor_len, 
    y = Inf, 
    label = paste0(" Overall median: ", format(round(median_anchor_len, 1), big.mark = ","), " bp"), 
    vjust = 2, 
    hjust = -0.05, 
    fontface = "bold", 
    color = "black"
  ) +
  labs(
    title = "Anchor length distribution (0 – 20,000 bp)",
    subtitle = paste0("Total number of anchors in range: ", format(nrow(dt_anchors_20k), big.mark = ","), 
                      " / ", format(nrow(dt_anchors), big.mark = ",")),
    x = "Anchor length (bp)",
    y = "Number of anchors in a bin"
  ) +
  custom_theme

pdf(output_histogram_anchors_20k, width = 7, height = 5)
print(p_hist_anchors_20k)
dev.off()

cat("\nSukces! Wygenerowano wszystkie wykresy i zapisano je w katalogu:\n", output_dir, "\n")
