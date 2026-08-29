library(data.table)
library(ggtern)
library(ggplot2)
library(viridis)
library(gridExtra)

##########################################################

n_sim <- 100000

rGDM3 <- function(n, fit, sorted_names) {
  alpha <- sapply(sorted_names[1:2], function(nm) fit@estimate[paste0("alpha_", nm)])
  beta  <- sapply(sorted_names[1:2], function(nm) fit@estimate[paste0("beta_", nm)])

  p1 <- rbeta(n, alpha[1], beta[1])
  p2 <- rbeta(n, alpha[2], beta[2]) * (1 - p1)
  p3 <- pmax(0, 1 - p1 - p2)

  res <- cbind(p1, p2, p3)
  colnames(res) <- sorted_names
  res
}

# Analityczna gęstość GDM3 w nieprzeskalowanej przestrzeni oryginalnej
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

safe_extract <- function(x, name) {
  if (is.null(x)) return(NA_real_)
  if (!(name %in% names(x))) return(NA_real_)
  as.numeric(x[name])
}

##########################################################

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop("Usage: Rscript plot_single_triangle_density.R output.pdf time_window loop_id [population]")
}

pdf_file <- args[1]
tw <- args[2]
loop_id <- args[3]

pop <- if (length(args) >= 4 && nzchar(args[4])) {
  as.integer(args[4])
} else {
  NA_integer_
}

cat("PDF:", pdf_file, "\n")
cat("Time window:", tw, "\n")
cat("Loop:", loop_id, "\n")
cat("Population:", ifelse(is.na(pop), "NOT PROVIDED", pop), "\n")

dir.create(dirname(pdf_file), recursive = TRUE, showWarnings = FALSE)

##########################################################

config_path <- file.path("config", "config.yml")
config <- yaml::yaml.load_file(config_path)

resultsdir <- config$paths$resultsdir

if (is.na(pop)) {
  input_dir_fit <- file.path(resultsdir, "MGLMfit_GDM", "real_data", "all", "init_1e-6")
  fit_file <- file.path(
    input_dir_fit,
    paste0("fit_real_", tw, "_", loop_id, ".rds")
  )
} else {
  input_dir_fit <- file.path(resultsdir, "MGLMfit_GDM", "real_data", "pops", "init_1e-6")
  fit_file <- file.path(
    input_dir_fit,
    paste0("fit_real_", tw, "_", loop_id, "_pop_", pop, ".rds")
  )
}

if (!file.exists(fit_file)) {
  stop("Fit file not found: ", fit_file)
}

##########################################################

seed_base <- 12345
set.seed(seed_base + ifelse(is.na(pop), 0, pop))

fit <- readRDS(fit_file)

if (is.null(fit)) {
  stop("Fit is NULL")
}

est <- fit@estimate
se <- fit@SE

alpha_pars <- names(est)[grep("^alpha_", names(est))]
fitted_names <- sub("^alpha_", "", alpha_pars)

if (length(fitted_names) < 2 || !("x_out" %in% fitted_names)) {
  stop("Invalid fitted names in model object")
}

other_fitted <- setdiff(fitted_names, "x_out")[1]
remaining_name <- setdiff(c("x_A1", "x_A2"), other_fitted)
sorted_names <- c("x_out", other_fitted, remaining_name)

# Symulacja pomocnicza do wyznaczenia średniej i wartości progu
y_sim <- rGDM3(n_sim, fit, sorted_names)
dt <- as.data.table(y_sim)

ct <- suppressWarnings(cor.test(dt$x_A1, dt$x_A2, method = "spearman"))
rho <- as.numeric(ct$estimate)

mean_out <- mean(dt$x_out)
threshold <- 2 * mean_out - 1

##########################################################
# GENEROWANIE SIATKI W PRZESTRZENI PRZESKALOWANEJ I PRZELICZENIE GĘSTOŚCI
##########################################################

grid_res <- 300
vals <- seq(0.0005, 0.9995, length.out = grid_res)
grid_dt <- as.data.table(expand.grid(x_A1 = vals, x_out = vals))
grid_dt[, x_A2 := 1 - x_A1 - x_out]
grid_dt <- grid_dt[x_A2 > 0]

# Odwrotne przeskalowanie współrzędnych do przestrzeni oryginalnej [0, 1]
grid_dt[, x_out_orig := x_out * (1 - threshold) + threshold] # x_out = (x_out_orig - threshold) / (1 - threshold) 
grid_dt[, x_A1_orig  := x_A1 * (1 - threshold)] # x_A1 = x_A1_orig / (1 - threshold)
grid_dt[, x_A2_orig  := x_A2 * (1 - threshold)] # x_A2 = x_A2_orig / (1 - threshold)

# Pobranie parametrów alfa i beta dla dopasowanego modelu
alpha_vec <- sapply(sorted_names[1:2], function(nm) fit@estimate[paste0("alpha_", nm)])
beta_vec  <- sapply(sorted_names[1:2], function(nm) fit@estimate[paste0("beta_", nm)])

# Mapowanie zmiennych na odpowiednie pozycje w modelu GDM3
val_1 <- grid_dt$x_out_orig
val_2 <- grid_dt[[paste0(other_fitted, "_orig")]]
val_3 <- grid_dt[[paste0(remaining_name, "_orig")]]

# Obliczenie analitycznej gęstości w punktach oryginalnych
grid_dt[, density := dgdm3(val_1, val_2, val_3, alpha_vec, beta_vec)]

##########################################################
# DIAGNOSTYKA I PRZYGOTOWANIE WSPÓŁRZĘDNYCH DLA GGTERN
##########################################################

cat("\n=== KONTROLA GĘSTOŚCI NA SIATCE ===\n")
cat("Liczba punktów w siatce:", nrow(grid_dt), "\n")
cat("Zakres wartości gęstości:", range(grid_dt$density, na.rm = TRUE), "\n")
cat("=======================================\n\n")

breaks_LR <- pretty(c(0, 1 - threshold))
breaks_LR <- breaks_LR[breaks_LR > 0 & breaks_LR <= 1 - threshold]
labels_LR <- sprintf("%g", breaks_LR)
labels_LR <- paste0("$", sub("-0", "-", sub("e(.*)", " \\\\cdot 10^{\\1}", labels_LR)), "$")

breaks_T <- pretty(c(threshold, 1))
breaks_T <- breaks_T[breaks_T >= threshold]
labels_T <- sub(" - 0$", "", sprintf("1 - %g", 1 - breaks_T))
labels_T <- paste0("$", sub("-0", "-", sub("e(.*)", " \\\\cdot 10^{\\1}", labels_T)), "$")

p1 <- ggtern(grid_dt, aes(x = x_A1, y = x_out, z = x_A2, color = density)) +
  geom_point(size = 0.5, stroke = 0) +
  labs(
    x = "p_A1",
    y = "p_out",
    z = "p_A2",
    subtitle = sprintf("Spearman's $\\rho = %.3f$", rho)
  ) +
  scale_L_continuous(breaks = breaks_LR / (1 - threshold), labels = labels_LR) +
  scale_T_continuous(breaks = (breaks_T - threshold) / (1 - threshold), labels = labels_T) +
  scale_R_continuous(breaks = breaks_LR / (1 - threshold), labels = labels_LR) +
  scale_color_viridis_c(
    option = "D",
    trans = "log2",
    na.value = "white",
    name = "Density"
  ) +
  theme_bw() +
  theme_showarrows() +
  theme_latex() +
  theme(tern.panel.expand = 0.6) +
  theme(tern.axis.arrow.sep = 0.4)

pdf(pdf_file, width = 5, height = 4)
print(p1)
dev.off()

cat("Saved density plot to:", pdf_file, "\n")
