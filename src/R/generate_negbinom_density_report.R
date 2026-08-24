library(data.table)
library(dplyr)
library(ggplot2)

config_path <- file.path("config", "config.yml")
config <- yaml::yaml.load_file(config_path)

resultsdir <- config$paths$resultsdir
total_reads_base_dir <- file.path(resultsdir, "counts", "total_reads")
negbinom_dir  <- file.path(resultsdir, "parameters", "negbinom_distribution_of_total_reads")
output_pdf <- file.path(negbinom_dir, "report_negbinom_densities_with_histograms.pdf")

time_windows <- c("00-02", "02-04", "04-06", "06-08", "08-10", 
                  "10-12", "12-14", "14-16", "16-18", "18-20")

# Set the X-axis range (total_reads per cell count) to plot the PMF curves
x_seq <- 0:10000

message("=== Starting PDF Report Generation ===")
message("Input directory: ", negbinom_dir)

pdf(output_pdf, width = 14, height = 8.5)

# -------------------------------------------------------------------------
# PAGE 1: GLOBAL OVERVIEW FOR ALL 10 TIME WINDOWS (POPULATION 'ALL')
# -------------------------------------------------------------------------
message("Generating Page 1 (Global overview of all time windows)...")
plot_data_p1 <- list()
raw_reads_all_tw <- list() # Store raw counts for Page 2 histogram

for (tw in time_windows) {
  file_path <- file.path(negbinom_dir, paste0("negbinom_parameters_", tw, "_all_and_by_population.tsv.gz"))
  raw_reads_path <- file.path(total_reads_base_dir, paste0("total_reads_", tw, ".tsv.gz"))
  
  if(!file.exists(file_path)) next

  dt <- fread(file_path)[population == "all"]
  if(nrow(dt) == 0 || is.na(dt$mu)) next
  
  mu <- dt$mu
  size <- dt$size
  
  # R's dnbinom uses (size, prob) parameterization. Convert from (mu, size):
  p <- size / (size + mu)
  y_vals <- dnbinom(x_seq, size = size, prob = p)
  
  lbl <- paste0(tw, " (mu=", round(mu, 1), ", size=", round(size, 2), ")")
  plot_data_p1[[tw]] <- data.frame(x = x_seq, y = y_vals, tw = tw, legend_lbl = lbl)
  
  # Load raw counts for the global histogram on page 2
  if (file.exists(raw_reads_path)) {
    raw_dt <- fread(raw_reads_path, select = c("population", "total_reads"))
    raw_dt[, tw := tw]
    raw_reads_all_tw[[tw]] <- raw_dt
  }
}

if(length(plot_data_p1) > 0) {
  df_p1 <- bind_rows(plot_data_p1)
  
  p1 <- ggplot(df_p1, aes(x = x, y = y, color = legend_lbl)) +
    geom_line(linewidth = 1.2, alpha = 0.85) +
    labs(
      title = "Negative Binomial Distributions Across All Time Windows",
      subtitle = "Comparison of the global population ('all') using Probability Mass Function (PMF) curves",
      x = "Total Reads per Cell",
      y = "Probability Density (PMF)",
      color = "Time Window (Parameters)"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 16, color = "#0f172a"),
      plot.subtitle = element_text(size = 11, color = "#475569"),
      legend.position = "right",
      legend.text = element_text(size = 9.5),
      panel.grid.minor = element_blank()
    ) +
    scale_color_brewer(palette = "Spectral")
  
  print(p1)
} else {
  message("WARNING: No data available to generate the first page overview!")
}

# -------------------------------------------------------------------------
# PAGE 2: HISTOGRAMS FOR ALL TIME WINDOWS (POPULATION 'ALL')
# -------------------------------------------------------------------------
if (length(raw_reads_all_tw) > 0) {
  message("Generating Page 2 (Global histograms across all time windows)...")
  df_raw_all <- bind_rows(raw_reads_all_tw)
  
  # Obliczamy liczbę komórek dla każdego okna czasowego
  tw_counts <- df_raw_all[, .(n_cells = .N), by = tw]
  tw_counts[, tw_label := paste0("TW: ", tw, " (N=", format(n_cells, big.mark=","), ")")]
  
  # Dołączamy czytelne etykiety do danych
  df_raw_all <- merge(df_raw_all, tw_counts, by = "tw")
  
  p2 <- ggplot(df_raw_all, aes(x = total_reads, fill = tw)) +
    geom_histogram(bins = 100, alpha = 0.75, color = "white", linewidth = 0.1) +
    facet_wrap(~ tw_label, scales = "free_y", ncol = 5) +
    labs(
      title = "Empirical Total Reads Histograms Across All Time Windows",
      subtitle = "Global cell population distributions with total cell counts (N)",
      x = "Total Reads per Cell",
      y = "Cell Count (Frequency)"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 16, color = "#0f172a"),
      plot.subtitle = element_text(size = 11, color = "#475569"),
      legend.position = "none",
      strip.background = element_rect(fill = "#f1f5f9", color = NA),
      strip.text = element_text(face = "bold", color = "#334155")
    )
  print(p2)
}

# -------------------------------------------------------------------------
# PAGES 3+: PER-TIME-WINDOW DECOMPOSITION (CURVES AND THEN HISTOGRAMS)
# -------------------------------------------------------------------------
for (tw in time_windows) {
  file_path <- file.path(negbinom_dir, paste0("negbinom_parameters_", tw, "_all_and_by_population.tsv.gz"))
  raw_reads_path <- file.path(total_reads_base_dir, paste0("total_reads_", tw, ".tsv.gz"))
  
  if(!file.exists(file_path)) {
    message("Skipped window ", tw, " - .tsv.gz parameter file not found")
    next
  }
  
  dt <- fread(file_path)
  message("Processing window ", tw, " (found populations: ", nrow(dt), ")...")
  
  # A) Theoretical Fit Curves Generation
  plot_data_tw <- list()
  for(i in seq_len(nrow(dt))) {
    row <- dt[i]
    if(is.na(row$mu) || is.na(row$size)) next
    
    p <- row$size / (row$size + row$mu)
    y_vals <- dnbinom(x_seq, size = row$size, prob = p)
    
    lbl <- paste0(row$population, " (mu=", round(row$mu, 1), ", size=", round(row$size, 2), ")")
    plot_data_tw[[row$population]] <- data.frame(x = x_seq, y = y_vals, pop = row$population, legend_lbl = lbl)
  }
  
  if(length(plot_data_tw) == 0) next
  
  df_tw <- bind_rows(plot_data_tw)
  df_tw$line_weight <- ifelse(df_tw$pop == "all", "all_group", "subgroup")
  
  p_tw_curves <- ggplot(df_tw, aes(x = x, y = y, color = legend_lbl, size = line_weight)) +
    geom_line(alpha = 0.85) +
    scale_size_manual(values = c("all_group" = 1.7, "subgroup" = 0.9), guide = "none") +
    labs(
      title = paste("Time Window:", tw, "— Theoretical Fits"),
      subtitle = "Comparison of fitted Negative Binomial distributions across cell populations",
      x = "Total Reads per Cell",
      y = "Probability Density (PMF)",
      color = "Population (Parameters)"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 16, color = "#0f172a"),
      plot.subtitle = element_text(size = 11, color = "#475569"),
      legend.position = "right",
      legend.text = element_text(size = 9),
      panel.grid.minor = element_blank()
    )
  print(p_tw_curves)
  
  # B) Empirical Histograms Generation for Current Time Window (Sorted & Counted)
  if (file.exists(raw_reads_path)) {
    raw_dt_tw <- fread(raw_reads_path, select = c("population", "total_reads"))
    
    # Tworzymy osobną kopię danych dla grupy 'all' przed modyfikacjami
    raw_dt_all <- copy(raw_dt_tw)[, population := "all"]
    df_hist_tw <- rbindlist(list(raw_dt_all, raw_dt_tw))
    
    # 1. Obliczamy liczebność (N) dla każdej populacji w tym oknie
    pop_sizes <- df_hist_tw[, .(n_cells = .N), by = population]
    
    # 2. Sortujemy populacje malejąco według liczby komórek
    # Chcemy, aby 'all' zawsze było na pierwszym miejscu, a reszta była posortowana malejąco
    pop_sizes[, is_all := ifelse(population == "all", 1, 0)]
    pop_sizes <- pop_sizes[order(-is_all, -n_cells)]
    
    # 3. Tworzymy nową etykietę łączącą nazwę i liczbę komórek (np. "1_AM (N=1,245)")
    pop_sizes[, pop_label := paste0(population, " (N=", format(n_cells, big.mark=","), ")")]
    
    # 4. Mapujemy nowe etykiety z powrotem do głównej tabeli i tworzymy czynnik (factor)
    # Użycie poziomu czynnika (levels) wymusza na ggplot2 dokładnie taką kolejność rysowania kafelków
    df_hist_tw <- merge(df_hist_tw, pop_sizes, by = "population")
    df_hist_tw[, pop_label := factor(pop_label, levels = pop_sizes$pop_label)]
    
    p_tw_hists <- ggplot(df_hist_tw, aes(x = total_reads, fill = population)) +
      geom_histogram(bins = 70, alpha = 0.8, color = "white", linewidth = 0.1) +
      facet_wrap(~ pop_label, scales = "free_y") +
      labs(
        title = paste("Time Window:", tw, "— Empirical Population Histograms"),
        subtitle = "Real counts distribution ordered by cell population size (N) from largest to smallest",
        x = "Total Reads per Cell",
        y = "Cell Count (Frequency)"
      ) +
      theme_minimal(base_size = 12) +
      theme(
        plot.title = element_text(face = "bold", size = 16, color = "#0f172a"),
        plot.subtitle = element_text(size = 11, color = "#475569"),
        legend.position = "none",
        strip.background = element_rect(fill = "#f8fafc", color = NA),
        strip.text = element_text(face = "bold", color = "#334155")
      )
    print(p_tw_hists)
  } else {
    message("WARNING: Missing raw data file for histogram profiling: ", raw_reads_path)
  }
}

# Close the PDF graphics device and save the file
dev.off()
message("=== Success: PDF Report generated successfully! ===")
message("Output file path: ", output_pdf)
