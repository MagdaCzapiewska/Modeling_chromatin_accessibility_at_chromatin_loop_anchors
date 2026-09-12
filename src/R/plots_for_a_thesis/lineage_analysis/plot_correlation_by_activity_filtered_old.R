library(data.table)
library(ggplot2)
library(patchwork)
library(stringr)
library(ggsignif)
library(ggtext)

# 1. Konfiguracja ścieżek
config <- yaml::yaml.load_file(file.path("config", "config.yml"))
datadir <- config$paths$datadir
resultsdir <- config$paths$resultsdir

loops_file <- file.path(datadir, "long_and_short_range_loops_D_mel.tsv")
loops_data <- fread(loops_file)

root_to_leaf_paths <- file.path(resultsdir, "lineage_analysis", "root_to_leaf_paths.txt")
output_file <- file.path(resultsdir, "plots_for_a_thesis", "lineage_analysis", "correlation_by_lineage_with_activity_filtered.pdf")
correlation_dir <- file.path(resultsdir, "MGLMfit_GDM_cor", "real_data", "pops", "init_1e-6")

time_windows <- c("00-02", "02-04", "04-06", "06-08", "08-10", "10-12", "12-14", "14-16", "16-18", "18-20")

target_columns <- c(
  "Dmel_6-8h_Neuroblasts", "Dmel_6-8h_Neurons", "Dmel_6-8h_Glia",
  "Dmel_10-12h_Neuroblasts", "Dmel_10-12h_Neurons", "Dmel_10-12h_Glia",
  "Dmel_14-16h_Neuroblasts", "Dmel_14-16h_Neurons", "Dmel_14-16h_Glia"
)

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

# 3. Filtrowanie wstępne i dołączanie danych z pętli
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

pdf(output_file, width = 16, height = 12)

for (p_idx in seq_along(paths)) {
  message("Przetwarzanie ścieżki: ", p_idx)
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
  
  # --- FILTROWANIE PĘTLI NA PODSTAWIE OKIEN RDZENIOWYCH (02-04 DO 16-18) ---
  core_tw <- c("02-04", "04-06", "06-08", "08-10", "10-12", "12-14", "14-16", "16-18")
  
  # Wyjątek dla ścieżki 18: ignorujemy 08-10 przy sprawdzaniu obecności
  if (p_idx == 18) {
    core_tw <- setdiff(core_tw, "08-10")
  }
  
  # Dynamiczny opis filtra w zależności od tego, czy jest to ścieżka 18
  filter_desc <- if (p_idx == 18) {
    "Fit quality (min EST/SE) >= 5 across core time windows (02-04 to 16-18, excluding 08-10)"
  } else {
    "Fit quality (min EST/SE) >= 5 across all core time windows (02-04 to 16-18)"
  }
  
  req_tw <- intersect(core_tw, path_dt$time_window)
  
  # Identyfikacja pętli obecnych we wszystkich wymaganych oknach rdzeniowych
  req_cor <- cor_path[time_window %in% req_tw & !is.na(loop_id) & !is.na(spearman_rho)]
  valid_loops <- req_cor[, .(n_tw = uniqueN(time_window)), by = loop_id][n_tw == length(req_tw), loop_id]
  
  # Zachowujemy te pętle dla WSZYSTKICH okien (w tym skrajnych 00-02 i 18-20)
  cor_path_filtered <- cor_path[loop_id %in% valid_loops]
  
  # Iterujemy przez 9 analizowanych kolumn
  for (t_col in target_columns) {
    
    plot_data <- copy(cor_path_filtered)
    plot_data[, split_var := factor(get(t_col), levels = c("0", "1"))]
    
    # Filtrowanie wartości NA z danych do wykresu, aby uniknąć problemów z osiowaniem punktów
    plot_data <- plot_data[!is.na(spearman_rho) & !is.na(split_var)]
    
    # Zliczanie N dla grup 0 i 1 po przefiltrowaniu
    counts_split <- plot_data[, .(n = .N), by = .(time_window, split_var)]
    counts_wide <- dcast(counts_split, time_window ~ split_var, value.var = "n", fill = 0, drop = FALSE)
    
    if (!"0" %in% names(counts_wide)) counts_wide[, `0` := 0]
    if (!"1" %in% names(counts_wide)) counts_wide[, `1` := 0]
    setnames(counts_wide, c("0", "1"), c("n_0", "n_1"))
    
    # Tworzenie etykiet osi X z naprzemiennym układem 2-poziomowym i formatowaniem liczności
    path_full <- merge(template_dt, path_dt, by = "time_window", all.x = TRUE)
    path_full <- merge(path_full, counts_wide, by = "time_window", all.x = TRUE)
    path_full[is.na(n_0), n_0 := 0]
    path_full[is.na(n_1), n_1 := 0]
    
    path_full <- path_full[order(x_base)]
    path_full[, pos_idx := .I]
    
    path_full[, x_label := ifelse(!is.na(population),
                                  ifelse(pos_idx %% 2 != 0,
                                         paste0(time_window, "\n", population, "\n\n", lineage, "\n(0: ", format_n(n_0), " | 1: ", format_n(n_1), ")"),
                                         paste0(time_window, "\n\n", population, "\n", lineage, "\n(0: ", format_n(n_0), " | 1: ", format_n(n_1), ")")),
                                  ifelse(pos_idx %% 2 != 0,
                                         paste0(time_window, "\n-\n\n-\n(0: 0 | 1: 0)"),
                                         paste0(time_window, "\n\n-\n-\n(0: 0 | 1: 0)")))]
    
    path_full$x_label <- factor(path_full$x_label, levels = path_full$x_label)
    
    plot_data <- merge(plot_data, path_full[, .(time_window, x_base, x_label)], by = "time_window", all.x = TRUE)
    
    clean_col_name <- gsub("Dmel_", "", t_col)
    clean_col_name <- gsub("_", " ", clean_col_name)
    
    # Funkcja tworząca połówkę wykresu wraz z testem Wilcoxona
    create_half_plot <- function(data_subset, labels_subset, title_text, hide_x_title = FALSE) {
      
      has_data <- nrow(data_subset[!is.na(spearman_rho) & !is.na(split_var)]) > 0
      
      p <- ggplot(data_subset, aes(x = x_label, y = spearman_rho, fill = split_var, color = split_var)) +
        scale_x_discrete(limits = labels_subset, drop = FALSE) +
        coord_cartesian(ylim = c(global_y_min, global_y_max + 0.35)) +
        scale_fill_manual(values = c("0" = "skyblue", "1" = "salmon"), drop = FALSE) +
        scale_color_manual(values = c("0" = "darkblue", "1" = "darkred"), drop = FALSE) +
        labs(
          title = title_text,
          x = if (hide_x_title) NULL else "Time window\nPopulation\nLineage\n(0: Inactive loops | 1: Active loops)",
          y = "Estimated Spearman's correlation (rho)",
          fill = paste("Activity status\n(", clean_col_name, ")"),
          color = paste("Activity status\n(", clean_col_name, ")")
        ) +
        theme_minimal(base_size = 12) +
        theme(
          plot.title = element_text(size = 13, face = "bold", color = "#2c3e50", margin = margin(t = 4, b = 4)),
          axis.title.x = element_text(size = 11.5, face = "bold", color = "#2c3e50", margin = margin(t = 10)),
          axis.title.y = element_text(size = 11.5, face = "bold", color = "#2c3e50", margin = margin(r = 8)),
          axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1, lineheight = 1.1, size = 10, color = "black"),
          axis.text.y = element_text(size = 10, color = "black"),
          panel.border = element_rect(color = "#dcdde1", fill = NA, linewidth = 0.4),
          panel.grid.minor = element_blank(),
          panel.grid.major.x = element_line(color = "#f1f2f6"),
          legend.title = element_text(size = 11, face = "bold", color = "#2c3e50"),
          legend.text = element_text(size = 10),
          legend.position = "right"
        )
      
      if (has_data) {
        p <- p + 
          geom_violin(
            position = position_dodge(width = 0.8),
            scale = "count", 
            alpha = 0.6, 
            color = "black", 
            na.rm = TRUE
          ) +
          geom_jitter(
            aes(group = split_var),
            position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8, seed = 42), 
            alpha = 0.4, 
            size = 0.8, 
            na.rm = TRUE
          )
        
        for (i in seq_along(labels_subset)) {
          lbl <- labels_subset[i]
          sub_dt <- data_subset[x_label == lbl & !is.na(spearman_rho) & !is.na(split_var)]
          
          n0 <- nrow(sub_dt[split_var == "0"])
          n1 <- nrow(sub_dt[split_var == "1"])
          
          if (n0 >= 3 && n1 >= 3) {
            wt <- wilcox.test(spearman_rho ~ split_var, data = sub_dt)
            p_val <- wt$p.value
            
            signif_stars <- ifelse(p_val < 0.001, "***",
                            ifelse(p_val < 0.01, "**",
                            ifelse(p_val < 0.05, "*", "ns")))
            
            label_text <- if (p_val < 0.001) "p < 0.001 (***)" else sprintf("p = %.3f (%s)", p_val, signif_stars)
            
            y_max_local <- max(sub_dt$spearman_rho, na.rm = TRUE)
            y_pos <- y_max_local + 0.08
            
            p <- p + geom_signif(
              annotation = label_text,
              xmin = i - 0.2,
              xmax = i + 0.2,
              y_position = y_pos,
              tip_length = 0.02,
              textsize = 3.8,
              color = "black"
            )
          }
        }
      } else {
        p <- p + annotate("text", x = 3, y = mean(c(global_y_min, global_y_max)), 
                          label = "No valid loops", color = "red", fontface = "italic", size = 4.5)
      }
      return(p)
    }
    
    labels_part1 <- path_full[x_base <= 5]$x_label
    labels_part2 <- path_full[x_base > 5]$x_label
    
    p1 <- create_half_plot(plot_data[x_base <= 5], labels_part1, 
                           paste("Time windows (00-02) - (08-10)"), hide_x_title = TRUE)
    
    p2 <- create_half_plot(plot_data[x_base > 5], labels_part2, 
                           paste("Time windows (10-12) - (18-20)"), hide_x_title = FALSE)
    
    final_plot <- p1 / p2 + 
      plot_annotation(
        title = bquote(bold("Estimated Spearman's correlation across cell lineage (Path " * .(p_idx) * " - Filtered)")),
        subtitle = paste0("<b>Split by:</b> Loop activity status (<i>", clean_col_name, "</i>) | ",
                          "<b>Difference measured by:</b> Wilcoxon rank-sum test<br>",
                          "<span style='color: #2c3e50; font-style: italic;'><b>Filter applied:</b> ", filter_desc, "</span>"),
        theme = theme(
          plot.title = element_text(size = 16, hjust = 0.5, color = "#2c3e50", margin = margin(t = 6, b = 4)),
          plot.subtitle = element_markdown(size = 11.5, hjust = 0.5, color = "#2c3e50", lineheight = 1.3, margin = margin(b = 6))
        )
      ) +
      plot_layout(guides = "collect")
    
    print(final_plot)
  }
}

dev.off()
cat("Sukces! Wygenerowano plik PDF ze wszystkimi przefiltrowanymi wykresami i testami Wilcoxona:\n", output_file, "\n")
