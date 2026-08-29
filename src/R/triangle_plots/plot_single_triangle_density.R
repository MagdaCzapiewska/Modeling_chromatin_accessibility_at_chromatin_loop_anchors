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

##########################################################

safe_extract <- function(x, name) {
  if (is.null(x)) return(NA_real_)
  if (!(name %in% names(x))) return(NA_real_)
  as.numeric(x[name])
}

##########################################################

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop("Usage: Rscript plot_single_triangle_counts.R output.pdf time_window loop_id [population]")
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

seed_base <- 12345
set.seed(seed_base + ifelse(is.na(pop), 0, pop))

fit <- readRDS(fit_file)

if (is.null(fit)) {
    stop("Fit is NULL")
}

# Wyznaczenie sorted_names na podstawie właściwości modelu
est <- fit@estimate
se  <- fit@SE

alpha_pars <- names(est)[grep("^alpha_", names(est))]
fitted_names <- sub("^alpha_", "", alpha_pars)

if (length(fitted_names) < 2 || !("x_out" %in% fitted_names)) {
  stop("Błąd: Niepoprawne nazwy zmiennych w obiekcie fit (brak x_out lub za mało parametrów alpha).")
}

other_fitted <- setdiff(fitted_names, "x_out")[1]
remaining_name <- setdiff(c("x_A1", "x_A2"), other_fitted)

sorted_names <- c("x_out", other_fitted, remaining_name)

# Korelacja Spearman'a do podtytułu
y_sim <- rGDM3(n_sim, fit, sorted_names)
dt_sim <- as.data.table(y_sim)
ct <- suppressWarnings(cor.test(dt_sim$x_A1, dt_sim$x_A2, method = "spearman"))
rho <- as.numeric(ct$estimate)
rm(y_sim, dt_sim)

##########################################################
# WYZNACZENIE ŚRODKÓW CIĘŻKOŚCI BINÓW (30 x 30)

N_bins <- 30

# 1. Trójkąty skierowane w górę (i + j + k = N_bins - 1)
up_list <- list()
for (i in 0:(N_bins - 1)) {
  for (j in 0:(N_bins - 1 - i)) {
    k <- N_bins - 1 - i - j
    p1 <- (3 * i + 1) / (3 * N_bins)
    p2 <- (3 * j + 1) / (3 * N_bins)
    p3 <- (3 * k + 1) / (3 * N_bins)
    up_list[[length(up_list) + 1]] <- data.table(p1 = p1, p2 = p2, p3 = p3)
  }
}

# 2. Trójkąty skierowane w dół (i + j + k = N_bins - 2)
down_list <- list()
if (N_bins >= 2) {
  for (i in 0:(N_bins - 2)) {
    for (j in 0:(N_bins - 2 - i)) {
      k <- N_bins - 2 - i - j
      p1 <- (3 * i + 2) / (3 * N_bins)
      p2 <- (3 * j + 2) / (3 * N_bins)
      p3 <- (3 * k + 2) / (3 * N_bins)
      down_list[[length(down_list) + 1]] <- data.table(p1 = p1, p2 = p2, p3 = p3)
    }
  }
}

grid_dt <- rbind(rbindlist(up_list), rbindlist(down_list))
setnames(grid_dt, c("p1", "p2", "p3"), sorted_names)

# Obliczenie gęstości w środkach binów
alpha1 <- safe_extract(fit@estimate, paste0("alpha_", sorted_names[1]))
beta1  <- safe_extract(fit@estimate, paste0("beta_", sorted_names[1]))
alpha2 <- safe_extract(fit@estimate, paste0("alpha_", sorted_names[2]))
beta2  <- safe_extract(fit@estimate, paste0("beta_", sorted_names[2]))

p1_v <- grid_dt[[sorted_names[1]]]
p2_v <- grid_dt[[sorted_names[2]]]
z2_v <- p2_v / (1 - p1_v)

eps <- 1e-6
z2_v <- pmin(pmax(z2_v, eps), 1 - eps)
p1_v_safe <- pmin(pmax(p1_v, eps), 1 - eps)

grid_dt[, log_density := dbeta(p1_v_safe, alpha1, beta1, log = TRUE) +
                         dbeta(z2_v, alpha2, beta2, log = TRUE) -
                         log(1 - p1_v_safe)]

grid_dt[, density := exp(log_density)]

##########################################################

mean_out <- alpha1 / (alpha1 + beta1)
threshold <- 2 * mean_out - 1

dts <- grid_dt[x_out >= threshold & is.finite(density) & density > 0, ]
dts[, x_out := (x_out - threshold) / (1 - threshold)]
dts[, x_A1 := x_A1 / (1 - threshold)]
dts[, x_A2 := x_A2 / (1 - threshold)]

# --- DIAGNOSTYKA ---
cat("\n=== KONTROLA GĘSTOŚCI W ŚRODKACH BINÓW (30 BINÓW) ===\n")
cat("Liczba punktów centralnych po progowaniu:", nrow(dts), "z", nrow(grid_dt), "\n")
cat("Zakres log-gęstości: [", min(dts$log_density), ",", max(dts$log_density), "]\n")
cat("======================================================\n\n")

breaks_LR <- pretty(c(0, 1 - threshold))
breaks_LR <- breaks_LR[breaks_LR > 0 & breaks_LR <= 1 - threshold]
labels_LR <- sprintf("%g", breaks_LR)
labels_LR <- paste0("$", sub("-0", "-", sub("e(.*)", " \\\\cdot 10^{\\1}", labels_LR)), "$")

breaks_T <- pretty(c(threshold, 1))
breaks_T <- breaks_T[breaks_T >= threshold]
labels_T <- sub(" - 0$", "", sprintf("1 - %g", 1 - breaks_T))
labels_T <- paste0("$", sub("-0", "-", sub("e(.*)", " \\\\cdot 10^{\\1}", labels_T)), "$")

p1 <- ggtern(dts, aes(x = x_A1, y = x_out, z = x_A2)) +
  geom_tri_tern(bins = N_bins, aes(weight = density, fill = after_stat(value))) +
  labs(
    x = "$p_{\\mathrm{A1}}$",
    y = "$p_{\\mathrm{out}}$",
    z = "$p_{\\mathrm{A2}}$",
    subtitle = sprintf("Spearman's $\\rho = %.3f$", rho),
    fill = "Density"
  ) +
  scale_L_continuous(breaks = breaks_LR / (1 - threshold), labels = labels_LR) +
  scale_T_continuous(breaks = (breaks_T - threshold) / (1 - threshold), labels = labels_T) +
  scale_R_continuous(breaks = breaks_LR / (1 - threshold), labels = labels_LR) +
  scale_fill_viridis_c(
    option = "D",
    trans = "log2",
    na.value = "white"
  ) +
  theme_bw() +
  theme_showarrows() +
  theme_latex() +
  theme(tern.panel.expand = 0.6) +
  theme(tern.axis.arrow.sep = 0.4)

pdf(pdf_file, width = 5, height = 4)
print(p1)
dev.off()

cat("Saved to:", pdf_file, "\n")
