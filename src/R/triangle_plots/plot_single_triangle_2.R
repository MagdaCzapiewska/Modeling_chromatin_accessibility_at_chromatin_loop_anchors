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
  stop("Usage: Rscript plot_single_triangle_2.R output.pdf time_window loop_id [population]")
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
pdf(pdf_file, width = 8, height = 8)

##########################################################

config_path <- file.path("config", "config.yml")
config <- yaml::yaml.load_file(config_path)

resultsdir <- config$paths$resultsdir
input_dir <- file.path(resultsdir, "dmel", "embrio", "p_fitting_random_loops", "all_loops_joint_all_pop_neural", "processed_data")

reads_file <- file.path(input_dir, paste0("reads_", tw, "_", loop_id, ".tsv"))
if (!file.exists(reads_file)) {
  stop("File missing:", reads_file)
}

reads <- fread(reads_file, select = c("x_A1", "x_A2", "total_reads"))
reads[, x_out := total_reads - x_A1 - x_A2]
n_cells <- nrow(reads)
  
col_sums <- colSums(reads[, .(x_A1, x_A2, x_out)])
ord <- order(col_sums, decreasing = TRUE)
sorted_names <- names(col_sums)[ord]
rm(reads)

if (is.na(pop)) {
  input_dir_fit <- file.path(resultsdir, "dmel", "embrio", "MGLM_GDM_fit_all")
  fit_file <- file.path(
    input_dir_fit,
    paste0("fit_", tw, "_", loop_id, "_all.rds")
  )
} else {
  input_dir_fit <- file.path(resultsdir, "dmel", "embrio", "MGLM_GDM_fit_pops")
  fit_file <- file.path(
    input_dir_fit,
    paste0("fit_", tw, "_", loop_id, "_pop_", pop, ".rds")
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
#threshold <- 2 * mean_out - 1
threshold <- 0.999

dt <- dt[x_out >= threshold]

#dt[, out_new := (x_out - threshold) / (1 - threshold)]

# A1 + A2 = 1 - out
# A1' + A2' = scale_factor * A1 + scale_factor * A2 = 1 - out'
# scale_factor = (1 - out') / (A1 + A2) = (1 - out') / (1 - out)
#dt[, scale_factor := (1 - out_new) / (1 - x_out)]

#dt[, `:=`(
#  x_A1 = x_A1 * scale_factor,
#  x_A2 = x_A2 * scale_factor,
#  x_out = out_new
#)]

#dt[, c("out_new", "scale_factor") := NULL]

##########################################################

#center[, out_new := (x_out - threshold) / (1 - threshold)]

#center[, scale_factor := (1 - out_new) / (1 - x_out)]

#center[, `:=`(
#  x_A1 = x_A1 * scale_factor,
#  x_A2 = x_A2 * scale_factor,
#  x_out = out_new
#)]

#center[, c("out_new", "scale_factor") := NULL]

##########################################################

p1 <- ggtern(dt, aes(x = x_A1, y = x_out, z = x_A2)) +
  geom_tri_tern(aes(fill = after_stat(count)),
                stat = "tri_tern", bins = 20000) +
  geom_point(data = center,
             aes(x = x_A1, y = x_out, z = x_A2),
             color = "magenta", size = 3) +
  scale_fill_viridis_c(
    option = "D",
    trans = "log2",
    na.value = "white"
  ) +
  theme_bw() +
  theme_showarrows() +
  tern_limits(T = 1.0, L = 1 - threshold, R = 1 - threshold) +
  labs(
    title = sprintf("Loop %s", loop_id),
    subtitle = sprintf(
      "%s | population %s\nSpearman = %.3f, threshold = %.6f",
      tw,
      ifelse(is.na(pop), "ALL", pop),
      rho,
      threshold
    )
  )

print(p1)

# https://github.com/tidyverse/ggplot2/issues/6374

#p2 <- ggtern(dt, aes(x = x_A1, y = x_out, z = x_A2)) +
  #stat_density_tern(
  #  aes(fill = after_stat(level)),
  #  geom = "polygon"
  #) +
#  geom_density_tern(aes(fill = after_stat(level)), bins = 20, h = c(1.0, 1.0), bdl = TRUE, bdl.val = 1e-20) +
#  geom_point(data = center,
#             aes(x = x_A1, y = x_out, z = x_A2),
#             color = "magenta", size = 3) +
#  scale_color_viridis_c(option = "A") +
#  theme_bw() +
#  theme_showarrows() +
#  labs(
#    title = sprintf("Loop %s (contour)", loop_id),
#    subtitle = sprintf(
#      "%s | population %s\nSpearman = %.3f, threshold = %.6f",
#      tw,
#      ifelse(is.na(pop), "ALL", pop),
#      rho,
#      threshold
#    )
#  )

#print(p2)
#grid.arrange(p1, p2, ncol = 1)

dev.off()

cat("Saved to:", pdf_file, "\n")
