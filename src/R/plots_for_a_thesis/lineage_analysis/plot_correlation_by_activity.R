library(data.table)
library(ggplot2)
library(patchwork)
library(stringr)
library(ggsignif)
library(ggtext)
library(yaml)

# 1. Konfiguracja ścieżek
config <- yaml::yaml.load_file(file.path("config", "config.yml"))
datadir <- config$paths$datadir
resultsdir <- config$paths$resultsdir

loops_file <- file.path(datadir, "long_and_short_range_loops_D_mel.tsv")
loops_data <- fread(loops_file)

root_to_leaf_paths <- file.path(resultsdir, "lineage_analysis", "root_to_leaf_paths.txt")
output_dir <- file.path(resultsdir, "plots_for_a_thesis", "lineage_analysis")
output_file <- file.path(output_dir, "correlation_by_lineage_with_activity.pdf")
stats_output_file <- file.path(output_dir, "correlation_by_activity_wilcoxon_results_bh_adjusted.tsv")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

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
cat("-> Wczytywanie plików korelacji...\n")
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
cat("-> Filtrowanie danych i dołączanie statusu aktywności pętli...\n")
all_cor <- all_cor[loop_id %in% loops_data$loop_id]
all_cor <- all_cor[min_est_over_se >= 5]

all_cor <- merge(all_cor, loops_data[, c("loop_id", target_columns), with = FALSE], 
                 by = "loop_id", all.x = TRUE)

global_y_min <- min(all_cor$spearman_rho, na.rm = TRUE)
global_y_max <- max(all_cor$spearman_rho, na.rm = TRUE)

# 4. Wczytanie ścieżek
paths <- readLines(root_to_leaf_paths)
paths <- paths[paths != ""]

template_dt <- data.table(time_window = time_windows, x_base = 1:10)

# ==============================================================================
# KROK 1: PASS 1 - TESTOWANIE UNIKALNYCH POPULACJI PER TARGET
# ==============================================================================
cat("\n=== KROK 1: Obliczanie testów Wilcoxona dla unikalnych populacji per target (Pass 1) ===\n")

# Wyciągnięcie unikalnych węzłów (time_window + population) ze wszystkich ścieżek
all_nodes_dt <- unique(rbindlist(lapply(paths, function(path_str) {
  path_nodes <- strsplit(path_str, "\\|")[[1]]
  path_nodes <- gsub("^\\(|\\)$", "", path_nodes)
  rbindlist(lapply(path_nodes, function(node) {
    parts <- strsplit(node, ";")[[1]]
    data.table(time_window = parts[1], population = parts[2], lineage = parts[3])
  }))
})), by = c("time_window", "population"))

test_results_list <- list()

for (t_col in target_columns) {
  for (i in seq_len(nrow(all_nodes_dt))) {
    tw <- all_nodes_dt$time_window[i]
    pop <- all_nodes_dt$population[i]
    lin <- all_nodes_dt$lineage[i]
    
    sub_dt <- all_cor[time_window == tw & population == pop]
    if (nrow(sub_dt) == 0) next
    
    sub_dt[, split_var := factor(get(t_col), levels = c("0", "1"))]
    sub_dt <- sub_dt[!is.na(spearman_rho) & !is.na(split_var)]
    
    n0 <- nrow(sub_dt[split_var == "0"])
    n1 <- nrow(sub_dt[split_var == "1"])
    
    if (n0 >= 3 && n1 >= 3 && length(unique(sub_dt$spearman_rho)) > 1) {
      wt <- tryCatch(wilcox.test(spearman_rho ~ split_var, data = sub_dt), error = function(e) NULL)
      if (!is.null(wt)) {
        test_results_list[[length(test_results_list) + 1]] <- data.table(
          target_col = t_col,
          time_window = tw,
          population = pop,
          lineage = lin,
          p_value = wt$p.value,
          n0 = n0,
          n1 = n1
        )
      }
    }
  }
}

p_val_results <- rbindlist(test_results_list)

# ==============================================================================
# KROK 2: GLOBALNA POPRAWKA BENJAMINI-HOCHBERGA (BH / FDR) I ZAPIS DO TSV
# ==============================================================================
cat("\n=== KROK 2: Zastosowanie poprawki Benjamini-Hochberga (BH) ===\n")
total_tests <- nrow(p_val_results)
cat(sprintf("-> Łączna liczba przetestowanych unikalnych kombinacji (M): %d\n", total_tests))

if (total_tests > 0) {
  p_val_results[, p_adj := p.adjust(p_value, method = "BH")]
  
  # Zapis pełnych unikalnych wyników statystycznych do pliku .tsv
  fwrite(p_val_results, stats_output_file, sep = "\t")
  cat(sprintf("-> Zapisano wyniki statystyczne do pliku TSV: %s\n", stats_output_file))
  
  sig_raw <- sum(p_val_results$p_value < 0.05, na.rm = TRUE)
  sig_adj <- sum(p_val_results$p_adj < 0.05, na.rm = TRUE)
  cat(sprintf("-> Liczba istotnych testów przed poprawką (p < 0.05): %d\n", sig_raw))
  cat(sprintf("-> Liczba istotnych testów po poprawce BH (p_adj < 0.05): %d\n", sig_adj))
} else {
  warning("Brak spełnionych warunków do wykonania testów Wilcoxona!")
}

# ==============================================================================
# KROK 3: PASS 2 - GENEROWANIE WYKRESÓW Z WYNIKAMI PO KOREKCJI BH
# ==============================================================================
cat("\n=== KROK 3: Generowanie wykresów i zapis do pliku PDF ===\n")

pdf(output_file, width = 16, height = 12)

for (p_idx in seq_along(paths)) {
  message(sprintf("[%d/%d] Generowanie wykresów dla Path %d...", p_idx, length(paths), p_idx))
  path_str <- paths[p_idx]
  
  path_nodes <- strsplit(path_str, "\\|")[[1]]
  path_nodes <- gsub("^\\(|\\)$", "", path_nodes)
  
  path_list <- lapply(path_nodes, function(node) {
    parts <- strsplit(node, ";")[[1]]
    data.table(time_window = parts[1], population = parts[2], lineage = parts[3])
  })
  path_dt <- rbindlist(path_list)
  cor_path <- merge(path_dt, all_cor, by = c("time_window", "population"), all.x = TRUE)
  
  for (t_col in target_columns) {
    
    plot_data <- copy(cor_path)
    plot_data[, split_var := factor(get(t_col), levels = c("0", "1"))]
    plot_data <- plot_data[!is.na(spearman_rho) & !is.na(split_var)]
    
    counts_split <- plot_data[, .(n = .N), by = .(time_window, split_var)]
    counts_wide <- dcast(counts_split, time_window ~ split_var, value.var = "n", fill = 0, drop = FALSE)
    
    if (!"0" %in% names(counts_wide)) counts_wide[, `0` := 0]
    if (!"1" %in% names(counts_wide)) counts_wide[, `1` := 0]
    setnames(counts_wide, c("0", "1"), c("n_0", "n_1"))
    
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
    
    # Funkcja generująca połowę wykresu
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
        theme_minimal(base_size = 13) +
        theme(
          plot.title         = element_text(size = 14, face = "bold", color = "black", margin = margin(t = 4, b = 4)),
          axis.title.x       = element_text(size = 13, face = "bold", color = "black", margin = margin(t = 10)),
          axis.title.y       = element_text(size = 13, face = "bold", color = "black", margin = margin(r = 8)),
          axis.text.x        = element_text(angle = 0, hjust = 0.5, vjust = 1, lineheight = 1.1, size = 11, color = "black"),
          axis.text.y        = element_text(size = 11, color = "black"),
          panel.border       = element_rect(color = "black", fill = NA, linewidth = 0.5),
          panel.grid.minor   = element_blank(),
          panel.grid.major.x = element_line(color = "#f1f2f6"),
          legend.title       = element_text(size = 12, face = "bold", color = "black"),
          legend.text        = element_text(size = 11, color = "black"),
          legend.position    = "right"
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
            size = 1.0, 
            na.rm = TRUE
          )
        
        # Nanoszenie klamer z p_adj po poprawce BH
        for (i in seq_along(labels_subset)) {
          lbl <- labels_subset[i]
          node_info <- path_full[x_label == lbl]
          tw_curr <- unique(node_info$time_window)
          pop_curr <- unique(node_info$population)
          
          if (length(tw_curr) == 1 && length(pop_curr) == 1 && !is.na(pop_curr)) {
            res_row <- p_val_results[target_col == t_col & time_window == tw_curr & population == pop_curr]
            
            if (nrow(res_row) == 1) {
              p_adj <- res_row$p_adj
              
              signif_stars <- ifelse(p_adj < 0.001, "***",
                              ifelse(p_adj < 0.01, "**",
                              ifelse(p_adj < 0.05, "*", "ns")))
              
              label_text <- if (p_adj < 0.001) "p_adj < 0.001 (***)" else sprintf("p_adj = %.3f (%s)", p_adj, signif_stars)
              
              sub_dt <- data_subset[x_label == lbl & !is.na(spearman_rho) & !is.na(split_var)]
              if (nrow(sub_dt) > 0) {
                y_max_local <- max(sub_dt$spearman_rho, na.rm = TRUE)
                y_pos <- y_max_local + 0.08
                
                p <- p + geom_signif(
                  annotation = label_text,
                  xmin = i - 0.2,
                  xmax = i + 0.2,
                  y_position = y_pos,
                  tip_length = 0.02,
                  textsize = 4.0,
                  color = "black"
                )
              }
            }
          }
        }
      } else {
        p <- p + annotate("text", x = 3, y = mean(c(global_y_min, global_y_max)), 
                          label = "No valid loops", color = "red", fontface = "italic", size = 5)
      }
      return(p)
    }
    
    labels_part1 <- path_full[x_base <= 5]$x_label
    labels_part2 <- path_full[x_base > 5]$x_label
    
    p1 <- create_half_plot(plot_data[x_base <= 5], labels_part1, 
                           "Time windows (00-02) - (08-10)", hide_x_title = TRUE)
    
    p2 <- create_half_plot(plot_data[x_base > 5], labels_part2, 
                           "Time windows (10-12) - (18-20)", hide_x_title = FALSE)
    
    final_plot <- p1 / p2 + 
      plot_annotation(
        title = bquote(bold("Estimated Spearman's correlation across cell lineage (Path " * .(p_idx) * ")")),
        subtitle = paste0("<b>Split by:</b> Loop activity status (<i>", clean_col_name, "</i>) | ",
                          "<b>Difference measured by:</b> Wilcoxon rank-sum test (BH FDR adjusted)<br>",
                          "<span style='color: black; font-style: italic;'><b>Filter applied:</b> Included chromatin loops with fit quality metric (min EST/SE) >= 5</span>"),
        theme = theme(
          plot.title    = element_text(size = 17, face = "bold", hjust = 0.5, color = "black", margin = margin(t = 6, b = 4)),
          plot.subtitle = element_markdown(size = 12.5, hjust = 0.5, color = "black", lineheight = 1.3, margin = margin(b = 6))
        )
      ) + 
      plot_layout(guides = "collect")
    
    print(final_plot)
  }
}

dev.off()
cat("\n=== SUKCES! Proces zakończony pomyślnie ===")
cat("\nWykresy ze skorygowanymi wartościami p_adj zapisano do pliku:\n", output_file)
cat("\nStatystyki i wartości p_adj zapisano do pliku:\n", stats_output_file, "\n")
