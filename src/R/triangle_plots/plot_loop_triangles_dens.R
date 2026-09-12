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
  stop("Użycie: Rscript plot_loop_triangles_dens.R <loop_id> <fit|reg>")
}

loop_id_arg <- as.character(args[1])
mode_type   <- tolower(args[2]) # "fit" lub "reg"

if (!mode_type %in% c("fit", "reg")) {
  stop("Parametr 'mode_type' musi przyjmować wartość 'fit' lub 'reg'.")
}

cat(sprintf("Uruchamianie generowania gęstości z plików eval: Loop = %s | Tryb = %s\n", loop_id_arg, mode_type))

##########################################################
# KONFIGURACJA I ŚCIEŻKI
##########################################################

config_path <- file.path("config", "config.yml")
if (!file.exists(config_path)) {
  stop("Nie znaleziono pliku konfiguracji: ", config_path)
}
config <- yaml::yaml.load_file(config_path)
resultsdir <- config$paths$resultsdir

# Katalog i plik docelowy PDF
out_dir <- file.path(resultsdir, "triangle_plots", "time", paste0("density_", mode_type))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
pdf_file <- file.path(out_dir, paste0("loop_", loop_id_arg, ".pdf"))

# Ścieżka do katalogu z wynikami ewaluacji okien czasowych
eval_dir <- file.path(resultsdir, "MGLMreg", "eval_0h+_time_default")

time_windows <- c("00-02", "02-04", "04-06", "06-08", "08-10", 
                  "10-12", "12-14", "14-16", "16-18", "18-20")

prefix <- paste0(mode_type, "_")
grid_res <- 300

##########################################################
# FUNKCJA ANALITYCZNEJ GĘSTOŚCI GDM3
##########################################################

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
  
  # Ekstrakcja parametrów rozkładu z pliku ewaluacji
  a_out <- row_dt[[paste0(prefix, "alpha_x_out_est")]]
  b_out <- row_dt[[paste0(prefix, "beta_x_out_est")]]
  
  a_A1  <- row_dt[[paste0(prefix, "alpha_x_A1_est")]]
  b_A1  <- row_dt[[paste0(prefix, "beta_x_A1_est")]]
  
  a_A2  <- row_dt[[paste0(prefix, "alpha_x_A2_est")]]
  b_A2  <- row_dt[[paste0(prefix, "beta_x_A2_est")]]
  
  if (is.na(a_out) || is.na(b_out) || a_out <= 0 || b_out <= 0) {
    cat(sprintf("Pominięto okno %s (brak poprawnych parametrów dla x_out)\n", tw))
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
  
  alpha_vec <- c(a_out, a2)
  beta_vec  <- c(b_out, b2)

  # Analityczna wartość oczekiwana p_out
  mean_out <- a_out / (a_out + b_out)
  threshold <- 2 * mean_out - 1
  #threshold <- threshold + (1 - threshold) * 999999999/10000000000

  # Generowanie siatki w przestrzeni przeskalowanej [0, 1]
  vals <- seq(0.0005, 0.9995, length.out = grid_res)
  grid_dt <- as.data.table(expand.grid(x_A1 = vals, x_out = vals))
  grid_dt[, x_A2 := 1 - x_A1 - x_out]
  grid_dt <- grid_dt[x_A2 > 0]

  # Odwrotne przeskalowanie współrzędnych do oryginalnej przestrzeni [0, 1]
  grid_dt[, x_out_orig := x_out * (1 - threshold) + threshold]
  grid_dt[, x_A1_orig  := x_A1 * (1 - threshold)]
  grid_dt[, x_A2_orig  := x_A2 * (1 - threshold)]

  # Wyznaczenie punktów dla gęstości analitycznej
  val_1 <- grid_dt$x_out_orig
  val_2 <- grid_dt[[paste0(cat2_name, "_orig")]]
  val_3 <- grid_dt[[paste0(cat3_name, "_orig")]]

  grid_dt[, density := dgdm3(val_1, val_2, val_3, alpha_vec, beta_vec)]

  # Podziałki i etykiety osi
  breaks_LR <- pretty(c(0, 1 - threshold))
  breaks_LR <- breaks_LR[breaks_LR > 0 & breaks_LR <= 1 - threshold]
  labels_LR <- sprintf("%g", breaks_LR)
  labels_LR <- paste0("$", sub("-0", "-", sub("e(.*)", " \\\\cdot 10^{\\1}", labels_LR)), "$")

  breaks_T  <- pretty(c(threshold, 1))
  breaks_T  <- breaks_T[breaks_T >= threshold]
  labels_T  <- sub(" - 0$", "", sprintf("1 - %g", 1 - breaks_T))
  labels_T  <- paste0("$", sub("-0", "-", sub("e(.*)", " \\\\cdot 10^{\\1}", labels_T)), "$")
  
  # Pobranie metryk z eval do nagłówków
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
  
  # Tworzenie wykresu gęstości
  p <- ggtern(grid_dt, aes(x = x_A1, y = x_out, z = x_A2, color = density)) +
    geom_point(size = 0.5, stroke = 0) +
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
    scale_color_viridis_c(option = "D", name = "Density") +
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
cat(sprintf("Zapisano wielostronicowy plik PDF z gęstością: %s\n", pdf_file))
