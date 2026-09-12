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

# population opcjonalne
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
  return(base_dt)
}

other_fitted <- setdiff(fitted_names, "x_out")[1]

remaining_name <- setdiff(c("x_A1", "x_A2"), other_fitted)

sorted_names <- c("x_out", other_fitted, remaining_name)

y_sim <- rGDM3(n_sim, fit, sorted_names)

dt <- as.data.table(y_sim)
ct <- suppressWarnings(cor.test(dt$x_A1, dt$x_A2, method = "spearman"))
rho <- as.numeric(ct$estimate)
pval <- ct$p.value
  
alpha_est <- setNames(rep(NA_real_,3), sorted_names)
beta_est  <- setNames(rep(NA_real_,3), sorted_names)
alpha_se  <- setNames(rep(NA_real_,3), sorted_names)
beta_se   <- setNames(rep(NA_real_,3), sorted_names)

for (i in seq_along(sorted_names)) {
  nm <- sorted_names[i]
  alpha_est[nm] <- safe_extract(fit@estimate, paste0("alpha_", nm))
  beta_est[nm]  <- safe_extract(fit@estimate, paste0("beta_", nm))
}

alpha_se[sorted_names[1]] <- fit@SE[[1]]
alpha_se[sorted_names[2]] <- fit@SE[[2]]
beta_se[sorted_names[1]]  <- fit@SE[[3]]
beta_se[sorted_names[2]]  <- fit@SE[[4]]

##########################################################

center <- dt[, .(
  x_A1 = mean(x_A1),
  x_A2 = mean(x_A2),
  x_out = mean(x_out)
)]

##########################################################

mean_out <- mean(dt$x_out)
# threshold = 1 - 2 * (1 - mean_out) = 2 * mean_out - 1
threshold <- 2 * mean_out - 1
#threshold <- threshold + (1 - threshold) * 999/1000

dts <- dt[x_out >= threshold, ]
dts[, x_out := (x_out - threshold) / (1 - threshold)]

# A1 + A2 = 1 - out
# A1' + A2' = scale_factor * A1 + scale_factor * A2 = 1 - out'
# scale_factor = (1 - out') / (A1 + A2) = (1 - out') / (1 - out)
# scale_factor = (1 - (out - threshold) / (1 - threshold)) / (1 - out) = 
# = ((1 - threshold - out + threshold) / (1 - threshold)) / (1 - out) = 
# = ((1 - out) / (1 - threshold)) / (1 - out) = 1 / (1 - threshold)
dts[, x_A1 := x_A1 / (1 - threshold)]
dts[, x_A2 := x_A2 / (1 - threshold)]

# --- DIAGNOSTYKA OSOBLIWOŚCI ---

cat("\n=== KONTROLA OSOBLIWOŚCI SYMULACJI ===\n")

# 1. Wartości niepoprawne (NA, NaN, Inf)
non_finite_count <- sum(!is.finite(as.matrix(dt)))
cat("Nienumerowane wartości (NA/NaN/Inf):", non_finite_count, "\n")

# 2. Osobliwości brzegowe (dokładne zera i jedynki)
zeros <- dt[, lapply(.SD, function(x) sum(x == 0))]
ones  <- dt[, lapply(.SD, function(x) sum(x == 1))]
cat("Liczba dokładnych zer (skupienie na brzegach):\n")
print(zeros)
cat("Liczba dokładnych jedynek:\n")
print(ones)

# 3. Poprawność sumy prawdopodobieństw (p1 + p2 + p3 == 1)
row_sums <- rowSums(dt)
sum_errors <- sum(abs(row_sums - 1) > 1e-6)
cat("Liczba wierszy, gdzie suma != 1 (błąd precyzji):", sum_errors, "\n")

# 4. Sprawdzenie parametrów Beta (czy gęstość dąży do nieskończoności przy 0/1)
cat("Parametry rozkładu Beta:\n")
for (nm in sorted_names[1:2]) {
  a <- fit@estimate[paste0("alpha_", nm)]
  b <- fit@estimate[paste0("beta_", nm)]
  if (!is.na(a) && a < 1) cat(sprintf("  UWAGA: alpha_%s = %.3f (< 1, gęstość ucieka do Inf przy 0)\n", nm, a))
  if (!is.na(b) && b < 1) cat(sprintf("  UWAGA: beta_%s = %.3f (< 1, gęstość ucieka do Inf przy 1)\n", nm, b))
}

# 5. Sprawdzenie filtrowania dts (czy dts nie jest pusty lub podzielony przez 0)
if (exists("dts")) {
  cat("\n=== KONTROLA FILTROVANIA GGTERN ===\n")
  cat("Współczynnik (1 - threshold):", 1 - threshold, "\n")
  cat("Liczba punktów po przeskalowaniu (nrow(dts)):", nrow(dts), " z ", nrow(dt), "\n")
  if (nrow(dts) == 0) warning("Ramka dts jest pusta! Progowanie wykluczyło wszystkie punkty.")
  if (abs(1 - threshold) < 1e-9) warning("Dzielenie przez zero przy skalowaniu dts!")
}
cat("=======================================\n\n")

breaks_LR <- pretty(c(0, 1 - threshold))
breaks_LR <- breaks_LR[breaks_LR > 0 & breaks_LR <= 1 - threshold]
labels_LR <- sprintf("%g", breaks_LR)
labels_LR <- paste0("$", sub("-0", "-", sub("e(.*)", " \\\\cdot 10^{\\1}", labels_LR)), "$")

breaks_T <- pretty(c(threshold, 1))
breaks_T <- breaks_T[breaks_T >= threshold]
labels_T <- sub(" - 0$", "", sprintf("1 - %g", 1 - breaks_T))
labels_T <- paste0("$", sub("-0", "-", sub("e(.*)", " \\\\cdot 10^{\\1}", labels_T)), "$")

p1 <- ggtern(dts, aes(x = x_A1, y = x_out, z = x_A2)) +
  geom_tri_tern(aes(fill = after_stat(count))) +
  # geom_point(data = center, ...) +
  labs(
    x = "p_A1$",
    y = "p_out",
    z = "p_A2",
    subtitle = sprintf("Spearman's $\\rho = %.3f$", rho)
  ) +
  scale_L_continuous(breaks = breaks_LR / (1 - threshold), labels = labels_LR) +
  scale_T_continuous(breaks = (breaks_T - threshold) / (1 - threshold), labels = labels_T) +
  scale_R_continuous(breaks = breaks_LR / (1 - threshold), labels = labels_LR) +
  scale_fill_viridis_c(
    option = "D"#,
    #trans = "log2",
    #na.value = "white"
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
