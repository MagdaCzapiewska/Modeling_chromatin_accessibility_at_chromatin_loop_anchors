library(data.table)
library(ggplot2)
library(stringr)
library(grid) # Potrzebne do obsługi jednostek unit() przy strzałce

# 1. Konfiguracja ścieżek
config <- yaml::yaml.load_file(file.path("config", "config.yml"))
datadir <- config$paths$datadir
resultsdir <- config$paths$resultsdir

loops_file <- file.path(datadir, "long_and_short_range_loops_D_mel.tsv")
loops_data <- fread(loops_file)

root_to_leaf_paths <- file.path(resultsdir, "lineage_analysis", "root_to_leaf_paths.txt")
output_dir <- file.path(resultsdir, "plots_for_a_thesis", "lineage_analysis")
output_file <- file.path(output_dir, "correlation_point_for_each_loop.pdf")

# Ścieżka do katalogu na osobne pliki i jego utworzenie
single_plots_dir <- file.path(output_dir, "correlation_point_for_each_loop")
dir.create(single_plots_dir, recursive = TRUE, showWarnings = FALSE)
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

# Ręczny JITTER (aby linie trafiały w kropki)
set.seed(42)
all_cor[, x_jittered := x_pos + runif(.N, min = -0.15, max = 0.15)]
all_cor[, y_jittered := spearman_rho + runif(.N, min = -0.005, max = 0.005)]

# Globalny zakres osi Y z dodatkowym zapasem na etykietę wyróżnionej pętli
global_y_min <- min(all_cor$y_jittered, na.rm = TRUE)
global_y_max <- max(all_cor$y_jittered, na.rm = TRUE)
y_margin <- (global_y_max - global_y_min) * 0.05
y_limits <- c(global_y_min - y_margin, global_y_max + y_margin + 0.06)

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
  
  # Wyciągnięcie danych TYLKO dla węzłów z tej ścieżki
  cor_path <- merge(path_dt, all_cor, by = c("time_window", "population"), all.x = TRUE)
  
  # Obliczenie ilości pętli w danych oknach
  counts <- cor_path[!is.na(spearman_rho), .(n_points = .N), by = time_window]
  
  # Tworzenie szkieletu osi X z uwzględnieniem pustych okien
  template_dt <- data.table(time_window = time_windows, x_base = 1:10)
  path_full <- merge(template_dt, path_dt, by = "time_window", all.x = TRUE)
  path_full <- merge(path_full, counts, by = "time_window", all.x = TRUE)
  path_full[is.na(n_points), n_points := 0]
  
  # Sortowanie po kolejności czasowej i indeksowanie pozycji
  path_full <- path_full[order(x_base)]
  path_full[, pos_idx := .I]
  
  # Tworzenie etykiet na oś X na 2 NAPRZEMIENNYCH POZIOMACH dla populacji
  path_full[, x_label := ifelse(!is.na(population),
                                ifelse(pos_idx %% 2 != 0,
                                       paste0(time_window, "\n", population, "\n\n", lineage, "\n(n=", n_points, ")"),
                                       paste0(time_window, "\n\n", population, "\n", lineage, "\n(n=", n_points, ")")),
                                ifelse(pos_idx %% 2 != 0,
                                       paste0(time_window, "\n-\n\n-\n(n=0)"),
                                       paste0(time_window, "\n\n-\n-\n(n=0)")))]
  
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
      title = paste0("Estimated Spearman's correlation across cell lineage (Path ", p_idx, ")"),
      subtitle = paste0("Points represent correlation for individual loops. Lines connect loop point to the point in the nearest valid ancestor in the path.\n",
                        "Filtered by fit quality (min_est_over_se) >= 5"),
      x = "Time window \n Population \n Lineage \n Number of loops (n)",
      y = "Estimated Spearman's correlation (rho)"
    ) +
    theme_bw(base_size = 13) +
    theme(
      plot.title         = element_text(size = 16, face = "bold", color = "black", hjust = 0.5),
      plot.subtitle      = element_text(size = 12, color = "black", lineheight = 1.25, hjust = 0.5, margin = margin(b = 8)),
      axis.title.x       = element_text(size = 13, face = "bold", color = "black", margin = margin(t = 10)),
      axis.title.y       = element_text(size = 13, face = "bold", color = "black", margin = margin(r = 8)),
      axis.text          = element_text(size = 10.5, color = "black"),
      axis.text.x        = element_text(angle = 0, hjust = 0.5, vjust = 1, lineheight = 1.1, size = 10.5, color = "black"),
      axis.text.y        = element_text(size = 10.5, color = "black"),
      panel.grid.minor   = element_blank(),
      panel.grid.major.x = element_line(color = "#f1f2f6")
    )
  
  if (has_data) {
    p <- p + 
      # Rysowanie linii łączących
      geom_line(data = plot_data, 
                aes(x = x_jittered, y = y_jittered, group = loop_id), 
                color = "gray60", alpha = 0.3, linewidth = 0.4) +
      # Rysowanie kropek
      geom_point(data = plot_data, 
                 aes(x = x_jittered, y = y_jittered), 
                 color = "darkblue", alpha = 0.6, size = 1.2)
    
    # WYRÓŻNIENIE PĘTLI O NAJWYŻSZEJ KORELACJI W OKNIE 08-10 DLA ŚCIEŻKI 19
    if (p_idx == 19) {
      dt_0810 <- plot_data[time_window == "08-10" & !is.na(spearman_rho)]
      if (nrow(dt_0810) > 0) {
        top_row <- dt_0810[which.max(spearman_rho)]
        top_loop_id <- top_row$loop_id
        top_x <- top_row$x_jittered
        top_y <- top_row$y_jittered
        
        p <- p +
          # Czerwone powiększone kółko na najwyższej kropce
          geom_point(data = top_row, 
                     aes(x = x_jittered, y = y_jittered), 
                     color = "red", size = 2.8) +
          # Czerwona strzałka skierowana w kropkę
          annotate("segment", 
                   x = top_x + 0.6, y = top_y + 0.07, 
                   xend = top_x + 0.06, yend = top_y + 0.015, 
                   color = "red", linewidth = 0.8, 
                   arrow = arrow(length = unit(0.25, "cm"), type = "closed")) +
          # Podpis z identyfikatorem loop_id
          annotate("text", 
                   x = top_x + 0.62, y = top_y + 0.075, 
                   label = paste0(top_loop_id), 
                   color = "red", fontface = "bold", size = 3.8, hjust = 0, vjust = 0)
      }
    }
  } else {
    p <- p + 
      annotate("text", x = 5.5, y = mean(y_limits), 
               label = "No loops with min_est_over_se >= 5", color = "red", size = 5)
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
