#!/usr/bin/env Rscript

library(data.table)
library(ggplot2)
library(ggtern)
library(viridis)
library(patchwork)
library(gridExtra)
library(grid)
library(stringr)

# ==============================================================================
# 1. Konfiguracja i obsługa argumentów wiersza poleceń
# ==============================================================================
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Użycie: Rscript plot_loop_time_analysis.R <loop_id> [output_pdf]")
}

loop_id_arg <- args[1]

config_path <- file.path("config", "config.yml")
if (!file.exists(config_path)) {
  if (file.exists("config.yml")) {
    config_path <- "config.yml"
  } else {
    stop("Nie znaleziono pliku konfiguracyjnego config.yml ani w folderze config/, ani w bieżącym katalogu.")
  }
}
config <- yaml::yaml.load_file(config_path)

resultsdir <- config$paths$resultsdir
srcdir     <- config$paths$srcdir

# Dedykowany katalog wyjściowy dla trójkątów w czasie
out_dir <- file.path(resultsdir, "triangle_plots", "time")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

output_pdf <- if (length(args) >= 2 && nzchar(args[2])) {
  args[2]
} else {
  file.path(out_dir, paste0("loop_", loop_id_arg, "_time_analysis.pdf"))
}

message("==================================================")
message(sprintf("Generowanie raportu z wyników eval dla pętli ID: %s", loop_id_arg))
message(sprintf("Plik docelowy: %s", output_pdf))
message("==================================================")

time_windows <- c("00-02", "02-04", "04-06", "06-08", "08-10", "10-12", "12-14", "14-16", "16-18", "18-20")
n_sim <- 50000

# Katalog zawierający pliki ewaluacyjne eval_cor_{tw}.tsv.gz
eval_dir <- file.path(resultsdir, "MGLMreg", "eval_0h+_time_default")

# ==============================================================================
# 2. Pomocnicze funkcje statystyczne i matematyczne
# ==============================================================================

# Symulacja Monte Carlo dla kaskadowego modelu GDM3
rGDM3_params <- function(n, sorted_names, alpha_vec, beta_vec) {
  if (is.null(alpha_vec) || is.null(beta_vec) || any(is.na(alpha_vec)) || any(is.na(beta_vec)) || any(alpha_vec <= 0) || any(beta_vec <= 0)) {
    return(NULL)
  }
  p1 <- rbeta(n, alpha_vec[1], beta_vec[1])
  p2 <- rbeta(n, alpha_vec[2], beta_vec[2]) * (1 - p1)
  p3 <- pmax(0, 1 - p1 - p2)
  
  res <- data.table(p1, p2, p3)
  setnames(res, sorted_names)
  return(res)
}

# Analityczna gęstość GDM3 w oryginalnej przestrzeni
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

# Ekstrakcja parametrów kaskady GDM z wiersza tabeli eval
extract_gdm_params_from_row <- function(row, prefix = "fit_") {
  alpha_out <- as.numeric(row[[paste0(prefix, "alpha_x_out_est")]])
  beta_out  <- as.numeric(row[[paste0(prefix, "beta_x_out_est")]])
  
  alpha_a1  <- as.numeric(row[[paste0(prefix, "alpha_x_A1_est")]])
  beta_a1   <- as.numeric(row[[paste0(prefix, "beta_x_A1_est")]])
  
  alpha_a2  <- as.numeric(row[[paste0(prefix, "alpha_x_A2_est")]])
  beta_a2   <- as.numeric(row[[paste0(prefix, "beta_x_A2_est")]])
  
  if (is.na(alpha_out) || is.na(beta_out) || alpha_out <= 0 || beta_out <= 0) return(NULL)
  
  if (!is.na(alpha_a1) && !is.na(beta_a1) && alpha_a1 > 0 && beta_a1 > 0) {
    step2_cat    <- "x_A1"
    leftover_cat <- "x_A2"
    alpha_step2  <- alpha_a1
    beta_step2   <- beta_a1
  } else if (!is.na(alpha_a2) && !is.na(beta_a2) && alpha_a2 > 0 && beta_a2 > 0) {
    step2_cat    <- "x_A2"
    leftover_cat <- "x_A1"
    alpha_step2  <- alpha_a2
    beta_step2   <- beta_a2
  } else {
    return(NULL)
  }
  
  sorted_names <- c("x_out", step2_cat, leftover_cat)
  alpha_vec    <- c(alpha_out, alpha_step2)
  beta_vec     <- c(beta_out,  beta_step2)
  
  list(
    sorted_names = sorted_names,
    alpha_vec    = alpha_vec,
    beta_vec     = beta_vec
  )
}

# ==============================================================================
# 3. Odczyt danych z wyników EVAL i przygotowanie prób okien czasowych
# ==============================================================================

time_summary_list <- list()
window_plots_data <- list()

for (tw in time_windows) {
  tw_bounds <- as.numeric(strsplit(tw, "-")[[1]])
  mid_time  <- mean(tw_bounds)
  
  eval_file <- file.path(eval_dir, paste0("eval_cor_", tw, ".tsv.gz"))
  
  if (!file.exists(eval_file)) {
    message(sprintf("Ostrzeżenie: Brak pliku eval dla okna %s (%s)", tw, eval_file))
    next
  }
  
  dt_eval <- fread(eval_file)
  row_eval <- dt_eval[loop_id == loop_id_arg]
  
  if (nrow(row_eval) == 0) {
    message(sprintf("Ostrzeżenie: Brak danych dla loop_id '%s' w oknie %s", loop_id_arg, tw))
    next
  }
  
  # --- 1. MGLMfit z danych eval ---
  fit_info <- extract_gdm_params_from_row(row_eval, prefix = "fit_")
  fit_dt   <- NULL
  fit_rho  <- as.numeric(row_eval$fit_spearman_rho)
  fit_means <- c(p_A1 = NA_real_, p_A2 = NA_real_, p_out = NA_real_)
  fit_alphas <- if (!is.null(fit_info)) fit_info$alpha_vec else NULL
  fit_betas  <- if (!is.null(fit_info)) fit_info$beta_vec else NULL
  fit_sorted_names <- if (!is.null(fit_info)) fit_info$sorted_names else c("x_out", "x_A1", "x_A2")
  
  if (!is.null(fit_info)) {
    set.seed(12345 + which(time_windows == tw))
    fit_dt <- rGDM3_params(n_sim, fit_sorted_names, fit_alphas, fit_betas)
    
    if (!is.null(fit_dt)) {
      if (is.na(fit_rho)) {
        ct <- suppressWarnings(cor.test(fit_dt$x_A1, fit_dt$x_A2, method = "spearman"))
        fit_rho <- as.numeric(ct$estimate)
      }
      fit_means <- c(p_A1 = mean(fit_dt$x_A1), p_A2 = mean(fit_dt$x_A2), p_out = mean(fit_dt$x_out))
    }
  }
  
  # --- 2. MGLMreg z danych eval ---
  reg_info <- extract_gdm_params_from_row(row_eval, prefix = "reg_")
  reg_dt   <- NULL
  reg_rho  <- as.numeric(row_eval$reg_spearman_rho)
  reg_means <- c(p_A1 = NA_real_, p_A2 = NA_real_, p_out = NA_real_)
  reg_alphas <- if (!is.null(reg_info)) reg_info$alpha_vec else NULL
  reg_betas  <- if (!is.null(reg_info)) reg_info$beta_vec else NULL
  reg_sorted_names <- if (!is.null(reg_info)) reg_info$sorted_names else c("x_out", "x_A1", "x_A2")
  
  if (!is.null(reg_info)) {
    set.seed(54321 + which(time_windows == tw))
    reg_dt <- rGDM3_params(n_sim, reg_sorted_names, reg_alphas, reg_betas)
    
    if (!is.null(reg_dt)) {
      if (is.na(reg_rho)) {
        ct <- suppressWarnings(cor.test(reg_dt$x_A1, reg_dt$x_A2, method = "spearman"))
        reg_rho <- as.numeric(ct$estimate)
      }
      reg_means <- c(p_A1 = mean(reg_dt$x_A1), p_A2 = mean(reg_dt$x_A2), p_out = mean(reg_dt$x_out))
    }
  }
  
  # Jeśli w pliku eval zapisano reg_pred_p_*, używamy preferowanych wartości z predict()
  pred_p_a1  <- as.numeric(row_eval$reg_pred_p_x_A1)
  pred_p_a2  <- as.numeric(row_eval$reg_pred_p_x_A2)
  pred_p_out <- as.numeric(row_eval$reg_pred_p_x_out)
  
  if (!is.na(pred_p_a1) && !is.na(pred_p_a2) && !is.na(pred_p_out)) {
    reg_means <- c(p_A1 = pred_p_a1, p_A2 = pred_p_a2, p_out = pred_p_out)
  }

  # Wspólny próg przycinania (threshold) dla zachowania spójności skali w oknie
  m_out_fit <- if (!is.null(fit_dt)) mean(fit_dt$x_out) else NA_real_
  m_out_reg <- if (!is.null(reg_dt)) mean(reg_dt$x_out) else NA_real_
  m_out_val <- max(c(m_out_fit, m_out_reg), na.rm = TRUE)
  common_threshold <- if (is.finite(m_out_val)) max(0, 2 * m_out_val - 1) else 0
  #common_threshold <- common_threshold + (1 - common_threshold) * 99/100

  time_summary_list[[tw]] <- data.table(
    time_window = tw,
    mid_time    = mid_time,
    fit_p_A1    = fit_means["p_A1"],
    fit_p_A2    = fit_means["p_A2"],
    fit_p_out   = fit_means["p_out"],
    fit_rho     = fit_rho,
    reg_p_A1    = reg_means["p_A1"],
    reg_p_A2    = reg_means["p_A2"],
    reg_p_out   = reg_means["p_out"],
    reg_rho     = reg_rho
  )

  window_plots_data[[tw]] <- list(
    tw = tw,
    mid_time = mid_time,
    threshold = common_threshold,
    fit_dt = fit_dt, fit_rho = fit_rho, fit_alphas = fit_alphas, fit_betas = fit_betas, fit_names = fit_sorted_names,
    reg_dt = reg_dt, reg_rho = reg_rho, reg_alphas = reg_alphas, reg_betas = reg_betas, reg_names = reg_sorted_names
  )
}

summary_dt <- rbindlist(time_summary_list)

if (nrow(summary_dt) == 0) {
  stop(sprintf("Brak danych ewaluacyjnych dla pętli ID: %s we wszystkich oknach czasowych.", loop_id_arg))
}

# ==============================================================================
# 4. STRONA 1: Zmiana prawdopodobieństw p_A1, p_A2, p_out oraz korelacji w czasie
# ==============================================================================

plot_trajectory <- function(dt, val_fit_col, val_reg_col, title_text, y_label) {
  dt_long <- melt(dt, id.vars = c("time_window", "mid_time"),
                  measure.vars = c(val_fit_col, val_reg_col),
                  variable.name = "Model", value.name = "Value")
  dt_long[, Model := factor(ifelse(Model == val_fit_col, "MGLMfit", "MGLMreg"), levels = c("MGLMfit", "MGLMreg"))]

  ggplot(dt_long, aes(x = mid_time, y = Value, color = Model, shape = Model)) +
    geom_line(linewidth = 1, alpha = 0.8) +
    geom_point(size = 3) +
    scale_x_continuous(breaks = summary_dt$mid_time, labels = summary_dt$time_window) +
    scale_color_manual(values = c("MGLMfit" = "#2b5c8f", "MGLMreg" = "#d95f02")) +
    labs(title = title_text, x = "Time Window", y = y_label) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
      plot.title  = element_text(size = 11, face = "bold"),
      legend.position = "top"
    )
}

p_traj_A1  <- plot_trajectory(summary_dt, "fit_p_A1",  "reg_p_A1",  "Probability Trajectory: p_A1",  "Probability (p_A1)")
p_traj_A2  <- plot_trajectory(summary_dt, "fit_p_A2",  "reg_p_A2",  "Probability Trajectory: p_A2",  "Probability (p_A2)")
p_traj_out <- plot_trajectory(summary_dt, "fit_p_out", "reg_p_out", "Probability Trajectory: p_out", "Probability (p_out)")
p_traj_rho <- plot_trajectory(summary_dt, "fit_rho",   "reg_rho",   "Spearman Correlation (rho)",     "Correlation (rho)")

page1_grid <- (p_traj_A1 | p_traj_A2) / (p_traj_out | p_traj_rho) +
  plot_annotation(
    title = paste("Loop ID:", loop_id_arg, "- Dynamics of GDM Parameters over Time (from eval)"),
    subtitle = "Comparison of discrete MGLMfit windows vs continuous MGLMreg trajectory",
    theme = theme(plot.title = element_text(size = 16, face = "bold"))
  )

# ==============================================================================
# 5. STRONY 2+: Histogramy MGLMfit vs MGLMreg na trójkątach (ggtern)
# ==============================================================================

make_ternary_histogram <- function(dt, sorted_names, rho_val, threshold, model_label, tw_name) {
  if (is.null(dt) || nrow(dt) == 0) {
    df_empty <- data.frame(x = 1/3, y = 1/3, z = 1/3)
    p <- ggtern(df_empty, aes(x = x, y = y, z = z)) +
      labs(title = paste0(tw_name, " | ", model_label), subtitle = "Brak danych modelu") +
      annotate("text", x = 1/3, y = 1/3, label = "BRAK DANYCH", color = "red", size = 5) +
      theme_bw() + theme_latex()
    return(p)
  }

  dts <- copy(dt[x_out >= threshold])
  if (threshold > 0 && threshold < 1) {
    dts[, x_out := (x_out - threshold) / (1 - threshold)]
    dts[, x_A1  := x_A1 / (1 - threshold)]
    dts[, x_A2  := x_A2 / (1 - threshold)]
  }

  breaks_LR <- pretty(c(0, 1 - threshold))
  breaks_LR <- breaks_LR[breaks_LR > 0 & breaks_LR <= 1 - threshold]
  labels_LR <- sprintf("%g", breaks_LR)
  labels_LR <- paste0("$", sub("-0", "-", sub("e(.*)", " \\\\cdot 10^{\\1}", labels_LR)), "$")

  breaks_T <- pretty(c(threshold, 1))
  breaks_T <- breaks_T[breaks_T >= threshold]
  labels_T <- sub(" - 0$", "", sprintf("1 - %g", 1 - breaks_T))
  labels_T <- paste0("$", sub("-0", "-", sub("e(.*)", " \\\\cdot 10^{\\1}", labels_T)), "$")

  rho_str <- if (is.na(rho_val)) "NA" else sprintf("%.3f", rho_val)

  p <- ggtern(dts, aes(x = x_A1, y = x_out, z = x_A2)) +
    geom_tri_tern(aes(fill = after_stat(count))) +
    labs(
      x = "$p_{A1}$", y = "$p_{out}$", z = "$p_{A2}$",
      title = paste0(tw_name, " | ", model_label),
      subtitle = sprintf("Spearman's $\\\\rho = %s$", rho_str)
    ) +
    scale_L_continuous(breaks = breaks_LR / (1 - threshold), labels = labels_LR) +
    scale_T_continuous(breaks = (breaks_T - threshold) / (1 - threshold), labels = labels_T) +
    scale_R_continuous(breaks = breaks_LR / (1 - threshold), labels = labels_LR) +
    #scale_fill_viridis_c(option = "D", trans = "log2", na.value = "white", name = "Count") +
    scale_fill_viridis_c(option = "D", name = "Count") +
    theme_bw() +
    theme_showarrows() +
    theme_latex() +
    theme(
      tern.panel.expand = 0.5,
      tern.axis.arrow.sep = 0.35,
      plot.title = element_text(size = 11, face = "bold"),
      plot.subtitle = element_text(size = 9.5, face = "italic")
    )
  return(p)
}

# ==============================================================================
# 6. STRONY N+: Analityczne gęstości GDM3 na trójkątach (ggtern)
# ==============================================================================

make_ternary_density <- function(sorted_names, alpha_vec, beta_vec, rho_val, threshold, model_label, tw_name) {
  if (is.null(alpha_vec) || is.null(beta_vec) || any(is.na(alpha_vec)) || any(is.na(beta_vec))) {
    df_empty <- data.frame(x = 1/3, y = 1/3, z = 1/3)
    p <- ggtern(df_empty, aes(x = x, y = y, z = z)) +
      labs(title = paste0(tw_name, " | ", model_label), subtitle = "Brak danych gęstości") +
      annotate("text", x = 1/3, y = 1/3, label = "BRAK DANYCH", color = "red", size = 5) +
      theme_bw() + theme_latex()
    return(p)
  }

  grid_res <- 100
  vals <- seq(0.001, 0.999, length.out = grid_res)
  grid_dt <- as.data.table(expand.grid(x_A1 = vals, x_out = vals))
  grid_dt[, x_A2 := 1 - x_A1 - x_out]
  grid_dt <- grid_dt[x_A2 > 0]

  # Odwrotne przekształcenie do oryginalnych współrzędnych
  grid_dt[, x_out_orig := x_out * (1 - threshold) + threshold]
  grid_dt[, x_A1_orig  := x_A1 * (1 - threshold)]
  grid_dt[, x_A2_orig  := x_A2 * (1 - threshold)]

  other_fitted   <- sorted_names[2]
  remaining_name <- sorted_names[3]

  val_1 <- grid_dt$x_out_orig
  val_2 <- grid_dt[[paste0(other_fitted, "_orig")]]
  val_3 <- grid_dt[[paste0(remaining_name, "_orig")]]

  grid_dt[, density := dgdm3(val_1, val_2, val_3, alpha_vec, beta_vec)]

  breaks_LR <- pretty(c(0, 1 - threshold))
  breaks_LR <- breaks_LR[breaks_LR > 0 & breaks_LR <= 1 - threshold]
  labels_LR <- sprintf("%g", breaks_LR)
  labels_LR <- paste0("$", sub("-0", "-", sub("e(.*)", " \\\\cdot 10^{\\1}", labels_LR)), "$")

  breaks_T <- pretty(c(threshold, 1))
  breaks_T <- breaks_T[breaks_T >= threshold]
  labels_T <- sub(" - 0$", "", sprintf("1 - %g", 1 - breaks_T))
  labels_T <- paste0("$", sub("-0", "-", sub("e(.*)", " \\\\cdot 10^{\\1}", labels_T)), "$")

  rho_str <- if (is.na(rho_val)) "NA" else sprintf("%.3f", rho_val)

  p <- ggtern(grid_dt, aes(x = x_A1, y = x_out, z = x_A2, color = density)) +
    geom_point(size = 1.0, stroke = 0) +
    labs(
      x = "$p_{A1}$", y = "$p_{out}$", z = "$p_{A2}$",
      title = paste0(tw_name, " | ", model_label),
      subtitle = sprintf("Spearman's $\\\\rho = %s$", rho_str)
    ) +
    scale_L_continuous(breaks = breaks_LR / (1 - threshold), labels = labels_LR) +
    scale_T_continuous(breaks = (breaks_T - threshold) / (1 - threshold), labels = labels_T) +
    scale_R_continuous(breaks = breaks_LR / (1 - threshold), labels = labels_LR) +
    #scale_color_viridis_c(option = "D", trans = "log2", na.value = "white", name = "Density") +
    scale_color_viridis_c(option = "D", name = "Density") +
    theme_bw() +
    theme_showarrows() +
    theme_latex() +
    theme(
      tern.panel.expand = 0.5,
      tern.axis.arrow.sep = 0.35,
      plot.title = element_text(size = 11, face = "bold"),
      plot.subtitle = element_text(size = 9.5, face = "italic")
    )
  return(p)
}

# ==============================================================================
# 7. GENEROWANIE RAPORTU PDF
# ==============================================================================

pdf(output_pdf, width = 14, height = 10)

# Strona 1: Trajektorie czasowe prawdopodobieństw i korelacji
print(page1_grid)

# Strony 2+: Histogramy na trójkątach (po 2 okna na stronę = 4 trójkąty)
for (page_i in seq(1, length(time_windows), by = 2)) {
  tws_sub <- time_windows[page_i:min(page_i + 1, length(time_windows))]
  
  hist_plots <- list()
  for (tw in tws_sub) {
    info <- window_plots_data[[tw]]
    if (is.null(info)) next
    
    p_fit <- make_ternary_histogram(info$fit_dt, info$fit_names, info$fit_rho, info$threshold, "MGLMfit", tw)
    p_reg <- make_ternary_histogram(info$reg_dt, info$reg_names, info$reg_rho, info$threshold, "MGLMreg", tw)
    
    hist_plots[[length(hist_plots) + 1]] <- p_fit
    hist_plots[[length(hist_plots) + 1]] <- p_reg
  }
  
  if (length(hist_plots) == 4) {
    grid.arrange(
      hist_plots[[1]], hist_plots[[2]],
      hist_plots[[3]], hist_plots[[4]],
      ncol = 2, nrow = 2,
      top = textGrob(paste("Loop ID:", loop_id_arg, "- Ternary Histograms (MGLMfit vs MGLMreg)"),
                     gp = gpar(fontsize = 15, font = 2))
    )
  } else if (length(hist_plots) >= 1) {
    grid.arrange(
      grobs = hist_plots,
      ncol = 2,
      top = textGrob(paste("Loop ID:", loop_id_arg, "- Ternary Histograms (MGLMfit vs MGLMreg)"),
                     gp = gpar(fontsize = 15, font = 2))
    )
  }
}

# Strony N+: Analityczne gęstości na trójkątach (po 2 okna na stronę = 4 trójkąty)
for (page_i in seq(1, length(time_windows), by = 2)) {
  tws_sub <- time_windows[page_i:min(page_i + 1, length(time_windows))]
  
  dens_plots <- list()
  for (tw in tws_sub) {
    info <- window_plots_data[[tw]]
    if (is.null(info)) next
    
    p_fit <- make_ternary_density(info$fit_names, info$fit_alphas, info$fit_betas, info$fit_rho, info$threshold, "MGLMfit", tw)
    p_reg <- make_ternary_density(info$reg_names, info$reg_alphas, info$reg_betas, info$reg_rho, info$threshold, "MGLMreg", tw)
    
    dens_plots[[length(dens_plots) + 1]] <- p_fit
    dens_plots[[length(dens_plots) + 1]] <- p_reg
  }
  
  if (length(dens_plots) == 4) {
    grid.arrange(
      dens_plots[[1]], dens_plots[[2]],
      dens_plots[[3]], dens_plots[[4]],
      ncol = 2, nrow = 2,
      top = textGrob(paste("Loop ID:", loop_id_arg, "- Analytical Densities (MGLMfit vs MGLMreg)"),
                     gp = gpar(fontsize = 15, font = 2))
    )
  } else if (length(dens_plots) >= 1) {
    grid.arrange(
      grobs = dens_plots,
      ncol = 2,
      top = textGrob(paste("Loop ID:", loop_id_arg, "- Analytical Densities (MGLMfit vs MGLMreg)"),
                     gp = gpar(fontsize = 15, font = 2))
    )
  }
}

dev.off()

cat("\n[Sukces] Wygenerowano raport PDF na podstawie wyników eval:\n", output_pdf, "\n")
