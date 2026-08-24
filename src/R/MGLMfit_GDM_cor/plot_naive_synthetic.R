library(data.table)
library(ggplot2)
library(patchwork)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript plot_naive_synthetic.R <input_tsv_gz> <expected_rho> <output_pdf>")
}

input_file   <- args[1]  
expected_rho <- as.numeric(args[2]) 
output_pdf   <- args[3]  

if (!file.exists(input_file)) stop(paste("Input file not found:", input_file))

dt <- fread(input_file)
if(!dir.exists(dirname(output_pdf))) dir.create(dirname(output_pdf), recursive = TRUE)

# Szeroki i wysoki layout dający optymalne proporcje dla mniejszej liczby wykresów na stronę
pdf(output_pdf, width = 22, height = 14)

make_syn_pair <- function(sub_dt, title_str, expected_val) {
  n_datasets <- nrow(sub_dt[!is.na(spearman_rho)])
  
  ps <- ggplot(sub_dt, aes(x = spearman_rho)) +
    geom_histogram(binwidth = 0.02, boundary = 0, fill = "#708090", color = "white", linewidth = 0.1) +
    xlim(-1, 1) +
    geom_vline(xintercept = expected_val, linetype = "dashed", color = "red", linewidth = 0.6) +
    labs(title = paste0(title_str, " [Spearman]"), subtitle = paste0("(datasets = ", n_datasets, ")"), x = "Rho", y = "Count") +
    theme_minimal() +
    theme(plot.title = element_text(size = 10, face = "bold"), plot.subtitle = element_text(size = 8), axis.title = element_text(size = 8))
  
  pp <- ggplot(sub_dt, aes(x = pearson_r)) +
    geom_histogram(binwidth = 0.02, boundary = 0, fill = "#4682B4", color = "white", linewidth = 0.1) +
    xlim(-1, 1) +
    geom_vline(xintercept = expected_val, linetype = "dashed", color = "red", linewidth = 0.6) +
    labs(title = paste0(title_str, " [Pearson]"), subtitle = paste0("(datasets = ", n_datasets, ")"), x = "R", y = "Count") +
    theme_minimal() +
    theme(plot.title = element_text(size = 10, face = "bold"), plot.subtitle = element_text(size = 8), axis.title = element_text(size = 8), axis.title.y = element_blank())
  
  return(ps + pp)
}

# ==============================================================================
# SEKCJA 1: Wykresy rozbite po N
# ==============================================================================
unique_n <- sort(unique(dt$n))
plots_n <- list()

for (current_n in unique_n) {
  sub_dt <- dt[n == current_n]
  plots_n[[as.character(current_n)]] <- make_syn_pair(sub_dt, paste0("n = ", current_n), expected_rho)
}

plot_chunks_n <- split(plots_n, ceiling(seq_along(plots_n) / 6))
for (i in seq_along(plot_chunks_n)) {
  combined_n <- wrap_plots(plot_chunks_n[[i]], ncol = 2, byrow = TRUE) + 
    plot_annotation(
      title = paste0("Naive Synthetic Data Correlation Histograms (by N) | Expected Rho = ", expected_rho),
      subtitle = paste0("Page ", i, " of ", length(plot_chunks_n), " | Red line: true rho. Left: Spearman, Right: Pearson."),
      theme = theme(plot.title = element_text(size = 14, face = "bold"))
    )
  print(combined_n)
}

# ==============================================================================
# SEKCJA 2: Wykresy rozbite po MU
# ==============================================================================
unique_mu <- sort(unique(dt$mu))
plots_mu <- list()

for (current_mu in unique_mu) {
  sub_dt <- dt[mu == current_mu]
  plots_mu[[as.character(current_mu)]] <- make_syn_pair(sub_dt, paste0("mu = ", current_mu), expected_rho)
}

combined_mu <- wrap_plots(plots_mu, ncol = 2, nrow = 3) + 
  plot_annotation(
    title = paste0("Naive Synthetic Data Correlations (by MU) | Expected Rho = ", expected_rho),
    subtitle = "Red dashed line: true rho. Left: Spearman, Right: Pearson.",
    theme = theme(plot.title = element_text(size = 14, face = "bold"))
  )
print(combined_mu)

# ==============================================================================
# SEKCJA 3: Wykresy rozbite po SIZE_NEGBINOM (Rozbicie na 2 strony po maks 6 obiektów)
# ==============================================================================
#unique_size <- unique(dt$size_nb)
#numeric_parts <- sort(as.numeric(unique_size[!unique_size %in% c("inf", "fixed")]))
#ordered_sizes <- c(unique_size[unique_size %in% c("inf", "fixed")], as.character(numeric_parts))
ordered_sizes <- c("inf", "fixed", "0.5", "1.0", "1.5", "2.0", "2.5", "3.0", "3.5", "10")

plots_size <- list()
for (current_size in ordered_sizes) {
  sub_dt <- dt[size_nb == current_size]
  plots_size[[current_size]] <- make_syn_pair(sub_dt, paste0("sizeNB = ", current_size), expected_rho)
}

# Dzielimy listę 12 modeli na porcje po max 6 modeli (czyli 12 pojedynczych wykresów na stronę w układzie 3x2)
plot_chunks_size <- split(plots_size, ceiling(seq_along(plots_size) / 6))

for (i in seq_along(plot_chunks_size)) {
  combined_size <- wrap_plots(plot_chunks_size[[i]], ncol = 2, nrow = 3, byrow = TRUE) + 
    plot_annotation(
      title = paste0("Naive Synthetic Data Correlations (by size_negbinom) | Expected Rho = ", expected_rho),
      subtitle = paste0("Size_nb Section - Page ", i, " of ", length(plot_chunks_size), " | Red line: true rho."),
      theme = theme(plot.title = element_text(size = 14, face = "bold"))
    )
  print(combined_size)
}

dev.off()
cat("Successfully generated spaced Naive Synthetic PDF report:", output_pdf, "\n")
