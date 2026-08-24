library(data.table)
library(ggplot2)
library(patchwork)
library(stringr)

# 1. Konfiguracja ścieżek
config <- yaml::yaml.load_file(file.path("config", "config.yml"))
datadir <- config$paths$datadir
resultsdir <- config$paths$resultsdir

loops_file <- file.path(datadir, "long_and_short_range_loops_D_mel.tsv")
loops_data <- fread(loops_file)

root_to_leaf_paths <- file.path(resultsdir, "lineage_analysis", "root_to_leaf_paths.txt")
output_file <- file.path(resultsdir, "lineage_analysis", "correlation_by_lineage_with_activity.pdf")
correlation_dir <- file.path(resultsdir, "MGLMfit_GDM_cor", "real_data", "pops", "init_1e-6")

time_windows <- c("00-02", "02-04", "04-06", "06-08", "08-10", "10-12", "12-14", "14-16", "16-18", "18-20")

target_columns <- c(
  "Dmel_6-8h_Neuroblasts", "Dmel_6-8h_Neurons", "Dmel_6-8h_Glia",
  "Dmel_10-12h_Neuroblasts", "Dmel_10-12h_Neurons", "Dmel_10-12h_Glia",
  "Dmel_14-16h_Neuroblasts", "Dmel_14-16h_Neurons", "Dmel_14-16h_Glia"
)

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

# 3. Filtrowanie i dołączanie danych z pętli
all_cor <- all_cor[loop_id %in% loops_data$loop_id]
all_cor <- all_cor[min_est_over_se >= 5]

all_cor <- merge(all_cor, loops_data[, c("loop_id", target_columns), with = FALSE], 
                 by = "loop_id", all.x = TRUE)

global_y_min <- min(all_cor$spearman_rho, na.rm = TRUE)
global_y_max <- max(all_cor$spearman_rho, na.rm = TRUE)

# 4. Wczytanie ścieżek
paths <- readLines(root_to_leaf_paths)
paths <- paths[paths != ""]

# Szablon osi X dla wygody
template_dt <- data.table(time_window = time_windows, x_base = 1:10)

pdf(output_file, width = 16, height = 12) # Zwiększono wysokość, by pomieścić dwa ułożone pionowo wykresy

for (p_idx in seq_along(paths)) {
  message(p_idx, "\n")
  path_str <- paths[p_idx]
  
  path_nodes <- strsplit(path_str, "\\|")[[1]]
  path_nodes <- gsub("^\\(|\\)$", "", path_nodes)
  
  path_list <- lapply(path_nodes, function(node) {
    parts <- strsplit(node, ";")[[1]]
    data.table(time_window = parts[1], population = parts[2], lineage = parts[3])
  })
  path_dt <- rbindlist(path_list)
  
  # Łączymy dane o przodkach z korelacjami
  cor_path <- merge(path_dt, all_cor, by = c("time_window", "population"), all.x = TRUE)
  
  # Iterujemy przez 9 analizowanych kolumn
  for (t_col in target_columns) {
    
    plot_data <- copy(cor_path)
    # Zabezpieczamy poziomy na wypadek braku jednej z grup
    plot_data[, split_var := factor(get(t_col), levels = c("0", "1"))]
    
    # --- NOWOŚĆ: Dokładne zliczanie N dla grup 0 i 1 oddzielnie ---
    counts_split <- plot_data[!is.na(spearman_rho), .(n = .N), by = .(time_window, split_var)]
    counts_wide <- dcast(counts_split, time_window ~ split_var, value.var = "n", fill = 0, drop = FALSE)
    
    # Zabezpieczenie na wypadek, gdyby dla danej ścieżki całkowicie brakowało np. grupy 0
    if (!"0" %in% names(counts_wide)) counts_wide[, `0` := 0]
    if (!"1" %in% names(counts_wide)) counts_wide[, `1` := 0]
    setnames(counts_wide, c("0", "1"), c("n_0", "n_1"))
    
    # Tworzenie kompletnych etykiet dla tej konkretnej zmiennej t_col
    path_full <- merge(template_dt, path_dt, by = "time_window", all.x = TRUE)
    path_full <- merge(path_full, counts_wide, by = "time_window", all.x = TRUE)
    path_full[is.na(n_0), n_0 := 0]
    path_full[is.na(n_1), n_1 := 0]
    
    path_full[, x_label := ifelse(!is.na(population),
                                  paste0(time_window, "\n", population, "\n", lineage, "\n(0: ", n_0, " | 1: ", n_1, ")"),
                                  paste0(time_window, "\n-\n-\n(0: 0 | 1: 0)"))]
    
    path_full <- path_full[order(x_base)]
    path_full$x_label <- factor(path_full$x_label, levels = path_full$x_label)
    
    # Łączymy z powrotem etykiety do danych wykresu
    plot_data <- merge(plot_data, path_full[, .(time_window, x_base, x_label)], by = "time_window", all.x = TRUE)
    
    clean_col_name <- gsub("Dmel_", "", t_col)
    clean_col_name <- gsub("_", " ", clean_col_name)
    
    # Funkcja generująca połowę wykresu, zapewniająca perfekcyjny dodge
    create_half_plot <- function(data_subset, labels_subset, title_text, hide_x_title = FALSE) {
      
      has_data <- nrow(data_subset[!is.na(spearman_rho) & !is.na(split_var)]) > 0
      
      # drop = FALSE w scale_x_discrete gwarantuje obecność wszystkich 5 kolumn
      p <- ggplot(data_subset, aes(x = x_label, y = spearman_rho, fill = split_var, color = split_var)) +
        scale_x_discrete(limits = labels_subset, drop = FALSE) +
        coord_cartesian(ylim = c(global_y_min, global_y_max)) +
        # drop = FALSE w paletach kolorów wymusza pozostawienie pustego miejsca, jeśli np. n_0 = 0
        scale_fill_manual(values = c("0" = "skyblue", "1" = "salmon"), drop = FALSE) +
        scale_color_manual(values = c("0" = "darkblue", "1" = "darkred"), drop = FALSE) +
        labs(
          title = title_text,
          x = if (hide_x_title) NULL else "Time window \n Population \n Lineage \n (0: Inactive loops | 1: Active loops)",
          y = "Spearman Rho",
          fill = paste("Activity status\n(", clean_col_name, ")"),
          color = paste("Activity status\n(", clean_col_name, ")")
        ) +
        theme_bw() +
        theme(
          axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1, lineheight = 1.1),
          legend.position = "right"
        )
      
      if (has_data) {
        # DODGING: kluczowe jest idealne zrównanie width = 0.8 w obu miejscach
        p <- p + 
          geom_violin(position = position_dodge(width = 0.8), scale = "count", 
                      alpha = 0.6, color = "black", na.rm = TRUE) +
          geom_jitter(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8, seed = 42), 
                      alpha = 0.4, size = 0.8, na.rm = TRUE)
      } else {
        p <- p + annotate("text", x = 3, y = mean(c(global_y_min, global_y_max)), 
                          label = "No valid loops", color = "red")
      }
      return(p)
    }
    
    # Dzielenie danych i etykiet na pierwszą (1:5) i drugą połowę (6:10) okien czasowych
    labels_part1 <- path_full[x_base <= 5]$x_label
    labels_part2 <- path_full[x_base > 5]$x_label
    
    p1 <- create_half_plot(plot_data[x_base <= 5], labels_part1, 
                           paste("Path", p_idx, "| Time windows 1-5"), hide_x_title = TRUE)
    
    p2 <- create_half_plot(plot_data[x_base > 5], labels_part2, 
                           paste("Path", p_idx, "| Time windows 6-10"), hide_x_title = FALSE)
    
    # --- NOWOŚĆ: Łączenie wykresów jeden pod drugim przy pomocy Patchwork ---
    final_plot <- p1 / p2 + 
      plot_annotation(
        title = paste("Complete Path:", p_idx),
        subtitle = paste("Split by loop activity in:", clean_col_name),
        theme = theme(plot.title = element_text(size = 14, face = "bold"))
      ) +
      plot_layout(guides = "collect") # scala wspólną legendę po prawej stronie
    
    print(final_plot)
  }
}

dev.off()
cat("Sukces! Zapisano podzielone wykresy do pliku:\n", output_file, "\n")
