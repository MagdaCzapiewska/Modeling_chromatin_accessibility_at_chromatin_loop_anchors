library(data.table)
library(ggplot2)
library(patchwork)
library(stringr)
library(ggtext)

# 1. Konfiguracja ścieżek
config <- yaml::yaml.load_file(file.path("config", "config.yml"))
datadir <- config$paths$datadir
resultsdir <- config$paths$resultsdir

loops_file <- file.path(datadir, "long_and_short_range_loops_D_mel.tsv")
loops_data <- fread(loops_file)

root_to_leaf_paths <- file.path(resultsdir, "lineage_analysis", "root_to_leaf_paths.txt")
output_dir <- file.path(resultsdir, "plots_for_a_thesis", "lineage_analysis")
output_file <- file.path(output_dir, "correlation_by_lineage_violin_filtered.pdf")

# Ścieżka do katalogu na osobne pliki i jego utworzenie
single_plots_dir <- file.path(output_dir, "correlation_by_lineage_violin_filtered")
dir.create(single_plots_dir, recursive = TRUE, showWarnings = FALSE)

correlation_dir <- file.path(resultsdir, "MGLMfit_GDM_cor", "real_data", "pops", "init_1e-6")

time_windows <- c("00-02", "02-04", "04-06", "06-08", "08-10", "10-12", "12-14", "14-16", "16-18", "18-20")

# Pomocnicza funkcja do formatowania liczności (np. 1,200)
format_n <- function(x) format(as.numeric(x), big.mark = ",", scientific = FALSE)

# 2. Wczytywanie danych korelacji
cor_list <- list()
for (tw in time_windows) {
  file_path <- file.path(correlation_dir, paste0("cor_", tw, ".tsv.gz"))
  if (file.exists(file_path)) {
    dt <- fread(file_path)
    cor_list[[tw]] <- dt
  } else {
    warning(paste("Plik nie istnieje:", file_path))
  }
}
all_cor <- rbindlist(cor_list)

# 3. Filtrowanie danych wg wytycznych
all_cor <- all_cor[loop_id %in% loops_data$loop_id]
all_cor <- all_cor[min_est_over_se >= 5]

# Obliczenie globalnego zakresu osi Y
global_y_min <- min(all_cor$spearman_rho, na.rm = TRUE)
global_y_max <- max(all_cor$spearman_rho, na.rm = TRUE)

# 4. Wczytanie ścieżek z pliku
paths <- readLines(root_to_leaf_paths)
paths <- paths[paths != ""]

# 5. Otwarcie zbiorczego pliku PDF do zapisu wszystkich wykresów
pdf(output_file, width = 14, height = 8)

for (p_idx in seq_along(paths)) {
  path_str <- paths[p_idx]
  
  # Rozbicie krotek wg separatora "|"
  path_nodes <- strsplit(path_str, "\\|")[[1]]
  path_nodes <- gsub("^\\(|\\)$", "", path_nodes)
  
  # Zbudowanie tabeli dla pojedynczej ścieżki
  path_list <- lapply(path_nodes, function(node) {
    parts <- strsplit(node, ";")[[1]]
    data.table(time_window = parts[1], population = parts[2], lineage = parts[3])
  })
  path_dt <- rbindlist(path_list)
  
  # Wyciągnięcie korelacji dla całej ścieżki
  cor_path <- merge(path_dt, all_cor, by = c("time_window", "population"), all.x = TRUE)
  
  # Wyznaczenie okien wymaganych i filtrowanie pętli
  req_tw <- c("02-04", "04-06", "06-08", "08-10", "10-12", "12-14", "14-16", "16-18")
  
  # Wyjątek dla Path 18: pominięcie okna 08-10 w sprawdzeniu
  if (p_idx == 18) {
    req_tw <- setdiff(req_tw, "08-10")
  }
  
  # Identyfikacja pętli, które występują we WSZYSTKICH wymaganych okienkach danej ścieżki
  req_cor <- cor_path[time_window %in% req_tw & !is.na(loop_id)]
  valid_loops <- req_cor[, .(n_tw = uniqueN(time_window)), by = loop_id][n_tw == length(req_tw), loop_id]
  
  # Przefiltrowanie cor_path tylko do wspólnych pętli
  cor_path <- cor_path[loop_id %in% valid_loops]
  
  # Liczenie punktów po przefiltrowaniu
  counts <- cor_path[!is.na(spearman_rho), .(n_points = .N), by = time_window]
  
  # Tworzenie szkieletu osi X
  template_dt <- data.table(time_window = time_windows)
  path_full <- merge(template_dt, path_dt, by = "time_window", all.x = TRUE)
  path_full <- merge(path_full, counts, by = "time_window", all.x = TRUE)
  path_full[is.na(n_points), n_points := 0]
  
  # Dodanie indeksu pozycji do układania populacji na 2 naprzemiennych wierszach
  path_full[, pos_idx := .I]
  
  path_full[, x_label := ifelse(!is.na(population),
                                ifelse(pos_idx %% 2 != 0,
                                       # Pozycja wyższa (dla nieparzystych): nazwa populacji, potem linia odstępu
                                       paste0(time_window, "\n", population, "\n\n", lineage, "\n(n=", format_n(n_points), ")"),
                                       # Pozycja niższa (dla parzystych): linia odstępu, potem nazwa populacji
                                       paste0(time_window, "\n\n", population, "\n", lineage, "\n(n=", format_n(n_points), ")")),
                                ifelse(pos_idx %% 2 != 0,
                                       paste0(time_window, "\n-\n\n-\n(n=0)"),
                                       paste0(time_window, "\n\n-\n-\n(n=0)")))]
  
  # Sztywne zablokowanie kolejności osi X
  path_full$x_label <- factor(path_full$x_label, levels = path_full$x_label)
  
  cor_plot <- merge(path_full[, .(time_window, x_label)], cor_path, by = "time_window", all.x = TRUE)
  
  # 6. Rysowanie wykresu
  has_data <- nrow(cor_plot[!is.na(spearman_rho)]) > 0
  
  p <- ggplot(cor_plot, aes(x = x_label, y = spearman_rho)) +
    scale_x_discrete(drop = FALSE) +
    coord_cartesian(ylim = c(global_y_min, global_y_max)) +
    labs(
      title = bquote(bold("Estimated Spearman's correlation along cell lineage - shared loops (Path " * .(p_idx) * ")")),
      subtitle = "<b>Filter applied:</b> Included chromatin loops with fit quality metric (min EST/SE ≥ 5) present across all required time windows",
      x = "Time window\nPopulation\nLineage\nNumber of loops (n)",
      y = "Estimated Spearman's correlation (rho)"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title = element_text(size = 13, hjust = 0.5, color = "#2c3e50", margin = margin(t = 6, b = 2)),
      plot.subtitle = element_markdown(size = 9.5, hjust = 0.5, color = "#2c3e50", lineheight = 1.25, margin = margin(b = 6)),
      axis.title.x = element_text(size = 10, face = "bold", color = "#2c3e50", margin = margin(t = 8)),
      axis.title.y = element_text(size = 10, face = "bold", color = "#2c3e50", margin = margin(r = 6)),
      axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1, lineheight = 1.1, size = 8.5, color = "black"),
      panel.border = element_rect(color = "#dcdde1", fill = NA, linewidth = 0.4),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_line(color = "#f1f2f6"),
      plot.margin = margin(t = 10, r = 10, b = 10, l = 10)
    )
  
  if (has_data) {
    p <- p + 
      geom_violin(scale = "count", fill = "skyblue", color = "black", alpha = 0.6, na.rm = TRUE) +
      geom_jitter(position = position_jitter(width = 0.1, seed = 42), 
                  alpha = 0.3, size = 0.8, color = "darkblue", na.rm = TRUE)
  } else {
    p <- p + 
      annotate("text", x = 5.5, y = mean(c(global_y_min, global_y_max)), 
               label = "No loops present across all required time windows with min EST/SE ≥ 5", color = "red", fontface = "italic")
  }
  
  # Zapis do zbiorczego pliku PDF
  print(p)
  
  # Zapis wykresu do osobnego pliku PDF
  single_pdf <- file.path(single_plots_dir, sprintf("path_%d.pdf", p_idx))
  ggsave(single_pdf, plot = p, width = 14, height = 8)
}

dev.off()

cat("Sukces!\n")
cat("1. Zbiorczy plik PDF:", output_file, "\n")
cat("2. Osobne pliki wygenerowane w katalogu:", single_plots_dir, "\n")
