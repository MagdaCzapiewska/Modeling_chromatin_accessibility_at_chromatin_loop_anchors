#!/usr/bin/env Rscript

library(data.table)
library(ggplot2)
library(patchwork)
library(stringr)
library(ggsignif)
library(ggtext)
library(yaml)

# ==============================================================================
# 1. Konfiguracja ścieżek i parametrów
# ==============================================================================
config <- yaml::yaml.load_file(file.path("config", "config.yml"))
datadir <- config$paths$datadir
resultsdir <- config$paths$resultsdir

loops_file <- file.path(datadir, "long_and_short_range_loops_D_mel.tsv")
activity_file <- file.path(resultsdir, "loop_activity_info.tsv.gz")

root_to_leaf_paths <- file.path(resultsdir, "lineage_analysis", "root_to_leaf_paths.txt")
output_dir <- file.path(resultsdir, "plots_for_a_thesis", "lineage_analysis")

correlation_dir <- file.path(resultsdir, "MGLMfit_GDM_cor", "real_data", "pops", "init_1e-6")

time_windows <- c("00-02", "02-04", "04-06", "06-08", "08-10", "10-12", "12-14", "14-16", "16-18", "18-20")

# ==============================================================================
# 2. Jednorazowe wczytanie danych korelacji i pętli
# ==============================================================================
message("Wczytywanie listy pętli oraz danych o korelacjach...")
loops_data <- fread(loops_file)

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
all_cor_base <- rbindlist(cor_list)

# Filtrowanie po dofitowaniu (min_est_over_se >= 5) oraz obecności w liście pętli
all_cor_base <- all_cor_base[loop_id %in% loops_data$loop_id & min_est_over_se >= 5]

paths <- readLines(root_to_leaf_paths)
paths <- paths[paths != ""]

# ==============================================================================
# 3. Definicje konfiguracji dla profili dokładnych oraz wzorcowych (star)
# ==============================================================================
configs <- list(
  exact = list(
    name = "exact",
    output_file = file.path(output_dir, "correlation_by_lineage_violin_profiles.pdf"),
    single_plots_dir = file.path(output_dir, "correlation_by_lineage_violin_profiles"),
    prof1_label = "Profile A (000010010)",
    prof2_label = "Profile B (000000010)",
    short_p1 = "A",
    short_p2 = "B",
    prof1_desc = "Profile A: present only in neurons in 10-12h and 14-16h",
    prof2_desc = "Profile B: present only in neurons in 14-16h",
    filter_func = function(dt) {
      dt[small_profile %in% c("000010010", "000000010"), 
         profile_group := ifelse(small_profile == "000010010", "Profile A (000010010)", "Profile B (000000010)")]
      return(dt[!is.na(profile_group)])
    }
  ),
  star = list(
    name = "star",
    output_file = file.path(output_dir, "correlation_by_lineage_violin_profiles_star.pdf"),
    single_plots_dir = file.path(output_dir, "correlation_by_lineage_violin_profiles_star"),
    prof1_label = "Profile A* (*0**1**1*)",
    prof2_label = "Profile B* (*0**0**1*)",
    short_p1 = "A*",
    short_p2 = "B*",
    prof1_desc = "Profile A*: present in neurons in 10-12h and 14-16h",
    prof2_desc = "Profile B*: present in neurons in 14-16h",
    filter_func = function(dt) {
      # *0**1**1* -> indeksy 2='0', 5='1', 8='1'
      # *0**0**1* -> indeksy 2='0', 5='0', 8='1'
      dt[grepl("^.0..1..1.$", small_profile), profile_group := "Profile A* (*0**1**1*)"]
      dt[grepl("^.0..0..1.$", small_profile), profile_group := "Profile B* (*0**0**1*)"]
      return(dt[!is.na(profile_group)])
    }
  )
)

# ==============================================================================
# 4. Funkcja generująca wykresy dla wybranej konfiguracji
# ==============================================================================
run_lineage_analysis <- function(cfg) {
  message(paste0("\n=================================================="))
  message(paste0(" Uruchamianie analizy dla trybu: ", cfg$name))
  message(paste0("=================================================="))
  
  dir.create(cfg$single_plots_dir, recursive = TRUE, showWarnings = FALSE)

  # Wczytanie i przygotowanie profilów aktywności
  act_dt <- fread(activity_file, select = c("loop_id", "small_profile"), colClasses = c(small_profile = "character"))
  act_dt[, small_profile := str_pad(as.character(small_profile), width = 9, side = "left", pad = "0")]

  # Przypisanie profilów wg reguły w konfiguracji
  act_dt <- cfg$filter_func(act_dt)
  act_dt[, profile_group := factor(profile_group, levels = c(cfg$prof1_label, cfg$prof2_label))]

  # Połączenie korelacji z przypisanymi profilami
  all_cor <- merge(all_cor_base, act_dt, by = "loop_id", all.x = FALSE)
  all_cor <- all_cor[!is.na(spearman_rho) & !is.na(profile_group)]

  if (nrow(all_cor) == 0) {
    warning("Brak danych po przefiltrowaniu dla trybu: ", cfg$name)
    return(NULL)
  }

  global_y_min <- min(all_cor$spearman_rho, na.rm = TRUE)
  global_y_max <- max(all_cor$spearman_rho, na.rm = TRUE)

  template_dt <- data.table(time_window = time_windows, x_base = 1:10)

  pdf(cfg$output_file, width = 15, height = 8)

  for (p_idx in seq_along(paths)) {
    message("  Przetwarzanie ścieżki ", p_idx, "/", length(paths))
    path_str <- paths[p_idx]
    
    path_nodes <- strsplit(path_str, "\\|")[[1]]
    path_nodes <- gsub("^\\(|\\)$", "", path_nodes)
    
    path_list <- lapply(path_nodes, function(node) {
      parts <- strsplit(node, ";")[[1]]
      data.table(time_window = parts[1], population = parts[2], lineage = parts[3])
    })
    path_dt <- rbindlist(path_list)
    
    cor_path <- merge(path_dt, all_cor, by = c("time_window", "population"), all.x = TRUE)
    cor_path <- cor_path[!is.na(spearman_rho) & !is.na(profile_group)]
    
    counts_split <- cor_path[, .(n = .N), by = .(time_window, profile_group)]
    counts_wide <- dcast(counts_split, time_window ~ profile_group, value.var = "n", fill = 0, drop = FALSE)
    
    if (!cfg$prof1_label %in% names(counts_wide)) counts_wide[, (cfg$prof1_label) := 0]
    if (!cfg$prof2_label %in% names(counts_wide)) counts_wide[, (cfg$prof2_label) := 0]
    setnames(counts_wide, c(cfg$prof1_label, cfg$prof2_label), c("n_p1", "n_p2"))
    
    path_full <- merge(template_dt, path_dt, by = "time_window", all.x = TRUE)
    path_full <- merge(path_full, counts_wide, by = "time_window", all.x = TRUE)
    path_full[is.na(n_p1), n_p1 := 0]
    path_full[is.na(n_p2), n_p2 := 0]
    
    path_full <- path_full[order(x_base)]
    path_full[, pos_idx := .I]
    
    # Etykiety osi X z uwzględnieniem odpowiednich skrótów profili (np. A/B lub A*/B*)
    path_full[, x_label := ifelse(!is.na(population),
                                  ifelse(pos_idx %% 2 != 0,
                                         paste0(time_window, "\n", population, "\n\n", lineage, 
                                                "\n(", cfg$short_p1, ": ", n_p1, " | ", cfg$short_p2, ": ", n_p2, ")"),
                                         paste0(time_window, "\n\n", population, "\n", lineage, 
                                                "\n(", cfg$short_p1, ": ", n_p1, " | ", cfg$short_p2, ": ", n_p2, ")")),
                                  ifelse(pos_idx %% 2 != 0,
                                         paste0(time_window, "\n-\n\n-\n(", cfg$short_p1, ": 0 | ", cfg$short_p2, ": 0)"),
                                         paste0(time_window, "\n\n-\n-\n(", cfg$short_p1, ": 0 | ", cfg$short_p2, ": 0)")))]
    
    path_full$x_label <- factor(path_full$x_label, levels = path_full$x_label)
    
    cor_plot <- merge(path_full[, .(time_window, x_base, x_label)], cor_path, by = "time_window", all.x = TRUE)
    cor_plot <- cor_plot[!is.na(spearman_rho) & !is.na(profile_group) & !is.na(x_label)]
    
    has_data <- nrow(cor_plot) > 0
    
    p <- ggplot(cor_plot, aes(x = x_label, y = spearman_rho, fill = profile_group, color = profile_group)) +
      scale_x_discrete(limits = levels(path_full$x_label), drop = FALSE) +
      coord_cartesian(ylim = c(global_y_min, global_y_max + 0.35)) +
      scale_fill_manual(values = setNames(c("#3498db", "#e74c3c"), c(cfg$prof1_label, cfg$prof2_label)), drop = FALSE) +
      scale_color_manual(values = setNames(c("#1f618d", "#922b21"), c(cfg$prof1_label, cfg$prof2_label)), drop = FALSE) +
      labs(
        title = paste0("Estimated Spearman's correlation across cell lineage (Path ", p_idx, ") - ", toupper(cfg$name)),
        subtitle = paste0("<b>Filtered by fit quality (min_est_over_se) &ge; 5</b><br>",
                          "<span style='font-size: 9.5pt; color: #444444;'>", cfg$prof1_desc, " &nbsp;|&nbsp; ", cfg$prof2_desc, "</span>"),
        x = "Time window\nPopulation\nLineage\n(Number of loops per profile)",
        y = "Spearman's rho",
        fill = "Profile Group",
        color = "Profile Group"
      ) +
      theme_bw(base_size = 11) +
      theme(
        plot.title = element_text(size = 14, face = "bold", color = "#2c3e50", hjust = 0.5),
        plot.subtitle = element_markdown(size = 10.5, color = "#2c3e50", hjust = 0.5, lineheight = 1.3),
        axis.title.x = element_text(size = 11, face = "bold", color = "#2c3e50", margin = margin(t = 8)),
        axis.title.y = element_text(size = 11, face = "bold", color = "#2c3e50", margin = margin(r = 8)),
        axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1, lineheight = 1.1, size = 9, color = "black"),
        axis.text.y = element_text(size = 9.5, color = "black"),
        legend.position = "top",
        legend.title = element_text(size = 10.5, face = "bold"),
        legend.text = element_text(size = 9.5),
        panel.grid.minor = element_blank()
      )
    
    if (has_data) {
      p <- p + 
        geom_violin(position = position_dodge(width = 0.8), scale = "count", 
                    alpha = 0.6, color = "black", na.rm = TRUE) +
        geom_jitter(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8, seed = 42), 
                    alpha = 0.4, size = 0.8, na.rm = TRUE)
      
      # Obliczanie testu Wilcoxona dla każdego okna czasowego
      labels_vec <- levels(path_full$x_label)
      for (i in seq_along(labels_vec)) {
        lbl <- labels_vec[i]
        sub_dt <- cor_plot[x_label == lbl & !is.na(spearman_rho) & !is.na(profile_group)]
        
        n1 <- nrow(sub_dt[profile_group == cfg$prof1_label])
        n2 <- nrow(sub_dt[profile_group == cfg$prof2_label])
        
        if (n1 >= 3 && n2 >= 3 && length(unique(sub_dt$spearman_rho)) > 1) {
          wt <- tryCatch(wilcox.test(spearman_rho ~ profile_group, data = sub_dt), error = function(e) NULL)
          
          if (!is.null(wt)) {
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
              textsize = 3.0,
              color = "black"
            )
          }
        }
      }
    } else {
      p <- p + annotate("text", x = 5.5, y = mean(c(global_y_min, global_y_max)), 
                        label = "No loops found for selected profiles", color = "red", fontface = "italic", size = 4.5)
    }
    
    print(p)
    
    single_pdf <- file.path(cfg$single_plots_dir, sprintf("path_%d.pdf", p_idx))
    ggsave(single_pdf, plot = p, width = 15, height = 8)
  }

  dev.off()
}

# ==============================================================================
# 5. Wykonanie analiz dla obu zestawów profili
# ==============================================================================
run_lineage_analysis(configs$exact)
run_lineage_analysis(configs$star)

cat("\nSukces! Wygenerowano wykresy dla profili dokładnych oraz star.\n")
