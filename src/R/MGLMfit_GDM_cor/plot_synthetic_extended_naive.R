library(data.table)
library(ggplot2)
library(patchwork)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript plot_synthetic_extended_naive.R <input_dir> <expected_rho> <output_base_dir>")
}

input_dir       <- args[1]
expected_rho    <- as.numeric(args[2])
output_base_dir <- args[3]

rho_str <- format(expected_rho, nsmall = 1)

pattern <- paste0("^naive_cor_seed.*\\.tsv\\.gz$")
files <- list.files(file.path(input_dir, paste0("rho_", rho_str)), pattern = pattern, full.names = TRUE)

if (length(files) == 0) {
  stop(paste("No correlation files found in:", file.path(input_dir, paste0("rho_", rho_str))))
}
dt <- rbindlist(lapply(files, fread))

# Pełne wektory extended
N_CELLS <- c("500", "1000", "1500", "2000", "2500", "3000", "3500", "4000", "5000", "10000", "20000", "40000", "60000", "80000", "100000", "120000", "140000", "160000")
MUS     <- c(1000, 2000, 3000, 4000, 5000, 6000)
SIZES   <- c("0.1", "0.2", "0.3", "0.4", "0.5", "0.6", "0.7", "0.8", "0.9", "1.0", "1.5", "2.0", "2.5", "3.0", "3.5", "10.0", "inf", "fixed")

dir.create(output_base_dir, recursive = TRUE, showWarnings = FALSE)

# Ulepszony helper: czytelne czcionki, marginesy i standardowe wielkości
plot_naive_pair <- function(sub_dt, title_txt) {
  n_datasets <- nrow(sub_dt[!is.na(spearman_rho)])
  
  ps <- ggplot(sub_dt, aes(x = spearman_rho)) +
    geom_histogram(binwidth = 0.04, boundary = 0, fill = "#708090", color = "white", linewidth = 0.1) +
    xlim(-1, 1) + 
    geom_vline(xintercept = expected_rho, linetype = "dashed", color = "red", linewidth = 0.6) +
    labs(title = paste0(title_txt, " [S]"), subtitle = paste0("ds=", n_datasets), x = "Rho", y = "Count") +
    theme_minimal() + 
    theme(
      plot.title = element_text(size = 10, face = "bold"), 
      plot.subtitle = element_text(size = 8), 
      axis.title = element_text(size = 8),
      panel.spacing = unit(0.5, "lines")
    )
  
  pp <- ggplot(sub_dt, aes(x = pearson_r)) +
    geom_histogram(binwidth = 0.04, boundary = 0, fill = "#4682B4", color = "white", linewidth = 0.1) +
    xlim(-1, 1) + 
    geom_vline(xintercept = expected_rho, linetype = "dashed", color = "red", linewidth = 0.6) +
    labs(title = paste0(title_txt, " [P]"), x = "R", y = "Count") +
    theme_minimal() + 
    theme(
      plot.title = element_text(size = 10, face = "bold"), 
      axis.title = element_text(size = 8), 
      axis.title.y = element_blank(),
      panel.spacing = unit(0.5, "lines")
    )
  
  return(ps + pp)
}

# --- PDF 1: Zmienia się N (Rozbicie 18 obiektów na 3 strony po 6 obiektów) ---
pdf(file.path(output_base_dir, "variable_N.pdf"), width = 20, height = 14)
for (current_mu in MUS) {
  for (current_size in SIZES) {
    plots <- list()
    for (current_n in N_CELLS) {
      sub_dt <- dt[mu == current_mu & size_nb == current_size & n == current_n]
      plots[[as.character(current_n)]] <- plot_naive_pair(sub_dt, paste0("n=", current_n))
    }
    
    # Podział 18 par na porcje po max 6 na stronę (układ 3 wiersze x 2 kolumny par)
    chunks <- split(plots, ceiling(seq_along(plots) / 6))
    for (page in seq_along(chunks)) {
      combined <- wrap_plots(chunks[[page]], ncol = 2, nrow = 3, byrow = TRUE) + 
        plot_annotation(
          title = paste0("Naive Extended | Variable N (mu=", current_mu, ", sizeNB=", current_size, ") | True rho = ", expected_rho),
          subtitle = paste0("Part ", page, " of ", length(chunks))
        )
      print(combined)
    }
  }
}
dev.off()

# --- PDF 2: Zmienia się MU (6 obiektów = idealnie 1 strona w układzie 3x2) ---
pdf(file.path(output_base_dir, "variable_MU.pdf"), width = 20, height = 14)
for (current_n in N_CELLS) {
  for (current_size in SIZES) {
    plots <- list()
    for (current_mu in MUS) {
      sub_dt <- dt[n == current_n & size_nb == current_size & mu == current_mu]
      plots[[as.character(current_mu)]] <- plot_naive_pair(sub_dt, paste0("mu=", current_mu))
    }
    
    # 6 obiektów mieści się idealnie na jednej stronie w siatce 3x2 par
    combined <- wrap_plots(plots, ncol = 2, nrow = 3, byrow = TRUE) + 
      plot_annotation(title = paste0("Naive Extended | Variable MU (n=", current_n, ", sizeNB=", current_size, ") | True rho = ", expected_rho))
    print(combined)
  }
}
dev.off()

# --- PDF 3: Zmienia się SIZE_NEGBINOM (Rozbicie 18 obiektów na 3 strony po 6 obiektów) ---
pdf(file.path(output_base_dir, "variable_SIZE.pdf"), width = 20, height = 14)
for (current_n in N_CELLS) {
  for (current_mu in MUS) {
    plots <- list()
    for (current_size in SIZES) {
      sub_dt <- dt[n == current_n & mu == current_mu & size_nb == current_size]
      plots[[current_size]] <- plot_naive_pair(sub_dt, paste0("size=", current_size))
    }
    
    # Podział 18 par na porcje po max 6 na stronę
    chunks <- split(plots, ceiling(seq_along(plots) / 6))
    for (page in seq_along(chunks)) {
      combined <- wrap_plots(chunks[[page]], ncol = 2, nrow = 3, byrow = TRUE) + 
        plot_annotation(
          title = paste0("Naive Extended | Variable SizeNB (n=", current_n, ", mu=", current_mu, ") | True rho = ", expected_rho),
          subtitle = paste0("Part ", page, " of ", length(chunks))
        )
      print(combined)
    }
  }
}
dev.off()

cat("Successfully generated spaced out Extended Naive PDF reports.\n")
