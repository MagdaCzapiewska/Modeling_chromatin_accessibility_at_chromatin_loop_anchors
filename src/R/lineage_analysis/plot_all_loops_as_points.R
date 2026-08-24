library(data.table)
library(ggplot2)
library(stringr)

# 1. Konfiguracja ścieżek
config <- yaml::yaml.load_file(file.path("config", "config.yml"))
datadir <- config$paths$datadir
resultsdir <- config$paths$resultsdir

loops_file <- file.path(datadir, "long_and_short_range_loops_D_mel.tsv")
loops_data <- fread(loops_file)

root_to_leaf_paths <- file.path(resultsdir, "lineage_analysis", "root_to_leaf_paths.txt")
output_file <- file.path(resultsdir, "lineage_analysis", "correlation_by_lineage_point_for_each_loop.pdf")
correlation_dir <- file.path(resultsdir, "MGLMfit_GDM_cor", "real_data", "pops", "init_1e-6")

time_windows <- c("00-02", "02-04", "04-06", "06-08", "08-10", "10-12", "12-14", "14-16", "16-18", "18-20")

# 2. Wczytywanie danych korelacji
cor_list <- list()
for (tw in time_windows) {
  file_path <- file.path(correlation_dir, paste0("cor_", tw, ".tsv.gz"))
  if (file.exists(file_path)) {
    cor_list[[tw]] <- fread(file_path)
  } else {
    warning(paste("Plik nie istnieje:", file_path))
  }
}
all_cor <- rbindlist(cor_list)

# 3. Filtrowanie i przygotowanie współrzędnych (Globalny Jitter)
all_cor <- all_cor[loop_id %in% loops_data$loop_id]
all_cor <- all_cor[min_est_over_se >= 5]

# Zamieniamy czas na pozycję numeryczną (1-10)
all_cor[, x_pos := match(time_window, time_windows)]

# --- Ręczny JITTER (aby linie trafiały w kropki) ---
set.seed(42)
all_cor[, x_jittered := x_pos + runif(.N, min = -0.15, max = 0.15)]
all_cor[, y_jittered := spearman_rho + runif(.N, min = -0.005, max = 0.005)]

# Globalny zakres osi Y
global_y_min <- min(all_cor$y_jittered, na.rm = TRUE)
global_y_max <- max(all_cor$y_jittered, na.rm = TRUE)
y_margin <- (global_y_max - global_y_min) * 0.05
y_limits <- c(global_y_min - y_margin, global_y_max + y_margin)

# 4. Wczytanie ścieżek z pliku
paths <- readLines(root_to_leaf_paths)
paths <- paths[paths != ""]

# 5. Otwarcie pliku PDF do zapisu wykresów
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
  
  # Wyciągnięcie danych TYLKO dla węzłów z tej ścieżki
  cor_path <- merge(path_dt, all_cor, by = c("time_window", "population"), all.x = TRUE)
  
  # Obliczenie ilości pętli w danych oknach
  counts <- cor_path[!is.na(spearman_rho), .(n_points = .N), by = time_window]
  
  # Tworzenie szkieletu osi X z uwzględnieniem pustych okien
  template_dt <- data.table(time_window = time_windows, x_base = 1:10)
  path_full <- merge(template_dt, path_dt, by = "time_window", all.x = TRUE)
  path_full <- merge(path_full, counts, by = "time_window", all.x = TRUE)
  path_full[is.na(n_points), n_points := 0]
  
  # Tworzenie etykiet na oś X
  path_full[, x_label := ifelse(!is.na(population),
                                paste0(time_window, "\n", population, "\n", lineage, "\n(n=", n_points, ")"),
                                paste0(time_window, "\n-\n-\n(n=0)"))]
  
  # Musimy posortować tabelę etykiet po czasie, by zachować kolejność na osi
  path_full <- path_full[order(x_base)]
  custom_x_labels <- path_full$x_label
  
  # Filtrujemy NA przed podaniem do ggplot, żeby linie dobrze się rysowały
  plot_data <- cor_path[!is.na(spearman_rho)]
  has_data <- nrow(plot_data) > 0
  
  # 6. Rysowanie wykresu
  p <- ggplot() +
    # Ustawiamy oś jako numeryczną ciągłą z miejscem na jitter (0.8 do 10.2)
    scale_x_continuous(breaks = 1:10, labels = custom_x_labels, limits = c(0.8, 10.2)) +
    coord_cartesian(ylim = y_limits) +
    labs(
      title = paste("Path", p_idx),
      subtitle = "Points represent individual loops. Lines connect the same loop to its nearest valid ancestor in the path.",
      x = "Time window \n Population \n Lineage \n Number of loops (n)",
      y = "Spearman Rho (min_est_over_se >= 5)"
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1, lineheight = 1.1),
      plot.subtitle = element_text(size = 10, color = "gray40")
    )
  
  if (has_data) {
    p <- p + 
      # Rysowanie linii łączących (group = loop_id sprawia, że każda pętla ma swoją linię!)
      geom_line(data = plot_data, 
                aes(x = x_jittered, y = y_jittered, group = loop_id), 
                color = "gray60", alpha = 0.3, linewidth = 0.4) +
      # Rysowanie kropek
      geom_point(data = plot_data, 
                 aes(x = x_jittered, y = y_jittered), 
                 color = "darkblue", alpha = 0.6, size = 1.2)
  } else {
    p <- p + 
      annotate("text", x = 5.5, y = mean(y_limits), 
               label = "No loops with min_est_over_se >= 5", color = "red")
  }
  
  print(p)
}

dev.off()

cat("Sukces! Wykresy (", length(paths), "stron ) zostały wygenerowane:\n", output_file, "\n")
