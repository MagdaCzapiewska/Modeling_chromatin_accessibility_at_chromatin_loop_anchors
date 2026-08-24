library(data.table)
library(ggplot2)
library(patchwork)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage: Rscript plot_real_correlations.R <mode> <init> <input_dir> <output_pdf> [<cardinality_file>]")
}

mode_arg   <- args[1]  # "all" lub "pops"
init_val   <- args[2]  # "default" lub "1e-6"
input_dir  <- args[3]  # katalog gdzie leżą pliki cor_{tw}.tsv.gz
output_pdf <- args[4]  # ścieżka do wynikowego pliku PDF

# Pobieramy ścieżkę kardynalności z 5. argumentu lub szukamy domyślnej

config_path <- file.path("config", "config.yml")

if (length(args) >= 5) {
  cardinality_file <- args[5]
} else if (file.exists(config_path)) {
  config <- yaml::yaml.load_file(config_path)
  resultsdir <- config$paths$resultsdir
  cardinality_file <- file.path(resultsdir, "cluster_cardinality.tsv.gz")
} else {
  # Awaryjny fallback, gdyby skrypt był odpalany całkiem poza Snakemake
  cardinality_file <- "results/cluster_cardinality.tsv.gz"
}

tw_list <- c("00-02", "02-04", "04-06", "06-08", "08-10", "10-12", "12-14", "14-16", "16-18", "18-20")

color_map <- c(
  "10" = "#99FF99",
  "5"  = "#FFFF99",
  "2"  = "#FFCC99",
  "1"  = "#FFB3B3",
  "0"  = "#99CCFF",
  "NA" = "#D3D3D3"
)

process_data <- function(dt) {
  param_cols <- c(
    "alpha_x_A1_est", "alpha_x_A2_est", "alpha_x_out_est",
    "beta_x_A1_est", "beta_x_A2_est", "beta_x_out_est",
    "alpha_x_A1_SE", "alpha_x_A2_SE", "alpha_x_out_SE",
    "beta_x_A1_SE", "beta_x_A2_SE", "beta_x_out_SE"
  )
  
  dt[, n_valid_params := rowSums(!is.na(.SD)), .SDcols = param_cols]
  
  dt[, color_group := cut(min_est_over_se, 
                          breaks = c(-Inf, 1, 2, 5, 10, Inf), 
                          labels = c("0", "1", "2", "5", "10"), 
                          right = FALSE)]
  
  dt[n_valid_params < 8 | is.na(min_est_over_se), color_group := "NA"]
  dt[, color_group := factor(color_group, levels = c("NA", "0", "1", "2", "5", "10"))]
  
  if ("population" %in% names(dt)) {
    dt[, pop_id := as.integer(sub("^([0-9]+)_.*", "\\1", population))]
    dt[, pop_name := population]
  }
  
  return(dt)
}

plot_type_1 <- function(dt, title_prefix = "", total_cells = NULL) {
  n_loops <- nrow(dt[!is.na(spearman_rho)])
  cell_info <- if(!is.null(total_cells)) paste0(", cells = ", total_cells) else ""
  full_title <- paste0(title_prefix, "\n(loops = ", n_loops, cell_info, ")")
  
  ggplot(dt, aes(x = spearman_rho, fill = color_group)) +
    geom_histogram(binwidth = 0.05, boundary = 0, color = "white", linewidth = 0.1) +
    scale_fill_manual(values = color_map, drop = FALSE, name = "min_est_over_se") +
    xlim(-1, 1) +
    labs(title = full_title, x = "Spearman Rho", y = "Count") +
    theme_minimal() +
    theme(plot.title = element_text(size = 9, face = "bold"), axis.title = element_text(size = 7))
}

plot_type_2 <- function(dt, title_prefix = "", total_cells = NULL) {
  dt_sub <- dt[as.numeric(as.character(color_group)) >= 5 & !is.na(as.numeric(as.character(color_group)))]
  n_loops <- nrow(dt_sub[!is.na(spearman_rho)])
  cell_info <- if(!is.null(total_cells)) paste0(", cells = ", total_cells) else ""
  full_title <- paste0(title_prefix, "\n(loops >= 5: ", n_loops, cell_info, ")")
  
  p <- ggplot(dt_sub, aes(x = spearman_rho)) +
    geom_density(fill = "#99FF99", alpha = 0.5) +
    xlim(-1, 1) +
    labs(title = full_title, subtitle = "min_est_over_se >= 5", x = "Spearman Rho", y = "Density") +
    theme_minimal() +
    theme(plot.title = element_text(size = 9, face = "bold"), axis.title = element_text(size = 7))
  
  if(nrow(dt_sub) < 2) {
    p = p + annotate("text", x = 0, y = 0.5, label = "No data (>=5)", size = 2)
  }
  return(p)
}

# --- Główny Nurt Wykonawczy ---

if(!dir.exists(dirname(output_pdf))) dir.create(dirname(output_pdf), recursive = TRUE)

# Wczytanie pliku kardynalności z miękkim lądowaniem (w razie braku pliku)
has_cardinality <- FALSE
if (file.exists(cardinality_file)) {
  dt_cardinality <- fread(cardinality_file)
  
  # Szukamy kolumny identyfikatora klastra (niezależnie od dokładnej pisowni)
  cluster_col <- grep("cluster", names(dt_cardinality), value = TRUE, ignore.case = TRUE)[1]
  
  if (length(cluster_col) > 0 && cluster_col %in% names(dt_cardinality)) {
    dt_card_long <- melt(dt_cardinality, 
                         id.vars = cluster_col, 
                         variable.name = "tw_window", 
                         value.name = "cells")
    
    setnames(dt_card_long, cluster_col, "cluster_id")
    dt_card_long[, cluster_id := as.integer(cluster_id)]
    dt_card_long[, tw_window := as.character(tw_window)]
    dt_card_long[, cells := as.integer(cells)]
    has_cardinality <- TRUE
  }
} else {
  warning(paste("Ostrzeżenie: Nie odnaleziono pliku liczności pod adresem:", cardinality_file, "- rysowanie bez liczby komórek."))
}

all_data_list <- list()
pdf(output_pdf, width = 16, height = 12)

for (tw in tw_list) {
  path <- file.path(input_dir, paste0("cor_", tw, ".tsv.gz"))
  if(!file.exists(path)) next
  
  dt <- fread(path)
  dt <- process_data(dt)
  dt[, tw_window := tw]
  all_data_list[[tw]] <- dt
  
  # Pobranie liczby komórek w zależności od trybu i dostępności pliku
  global_cells <- NULL
  if (mode_arg == "all" && has_cardinality) {
    global_cells <- sum(dt_card_long[tw_window == tw, cells], na.rm = TRUE)
  }
  
  if (mode_arg == "all") {
    p1 <- plot_type_1(dt, paste("TW:", tw), total_cells = global_cells)
    p2 <- plot_type_2(dt, paste("TW:", tw), total_cells = global_cells)
    print(p1 + p2)
    
  } else {
    # --- Tryb populacyjny ---
    pop_info <- unique(dt[!is.na(population), .(pop_id, pop_name)])
    
    if (has_cardinality) {
      pop_info <- merge(pop_info, 
                        dt_card_long[tw_window == tw, .(pop_id = cluster_id, cells)], 
                        by = "pop_id", 
                        all.x = TRUE)
    } else {
      pop_info[, cells := as.integer(NA)]
    }
    
    setorder(pop_info, pop_id) 
    if(nrow(pop_info) == 0) next
    
    p1_list <- list()
    p2_list <- list()
    
    for(i in seq_len(nrow(pop_info))) {
      pid   <- pop_info$pop_id[i]
      pname <- pop_info$pop_name[i]
      pcell <- if(has_cardinality) pop_info$cells[i] else NULL
      sub_dt <- dt[pop_id == pid]
      
      p1_list[[i]] <- plot_type_1(sub_dt, pname, total_cells = pcell) + theme(legend.position = "none")
      p2_list[[i]] <- plot_type_2(sub_dt, pname, total_cells = pcell)
    }
    
    pops_per_page <- 10
    plot_chunks_p1 <- split(p1_list, ceiling(seq_along(p1_list) / pops_per_page))
    plot_chunks_p2 <- split(p2_list, ceiling(seq_along(p2_list) / pops_per_page))
    
    for (page_idx in seq_along(plot_chunks_p1)) {
      print(wrap_plots(plot_chunks_p1[[page_idx]], ncol = 5) + 
            plot_layout(guides = "collect") +
            plot_annotation(
              title = paste("TW:", tw, "- Histograms per Population (Sorted)"),
              subtitle = paste0("Page ", page_idx, " of ", length(plot_chunks_p1))
            ))
    }
    
    for (page_idx in seq_along(plot_chunks_p2)) {
      print(wrap_plots(plot_chunks_p2[[page_idx]], ncol = 5) + 
            plot_layout(guides = "collect") +
            plot_annotation(
              title = paste("TW:", tw, "- Density plots (min_est_over_se >= 5)"),
              subtitle = paste0("Page ", page_idx, " of ", length(plot_chunks_p2))
            ))
    }
  }
}

if (mode_arg == "all" && length(all_data_list) > 0) {
  big_dt <- rbindlist(all_data_list)
  
  print(ggplot(big_dt, aes(x = spearman_rho, fill = color_group)) +
          geom_histogram(binwidth = 0.05, boundary = 0, color = "white", linewidth = 0.1) +
          scale_fill_manual(values = color_map, drop = FALSE, name = "min_est_over_se") +
          facet_wrap(~tw_window, nrow = 2) +
          labs(title = "Summary of All Time Windows - Histograms") + 
          theme_minimal())
  
  big_dt_sub <- big_dt[as.numeric(as.character(color_group)) >= 5]
  print(ggplot(big_dt_sub, aes(x = spearman_rho)) + 
          geom_density(fill = "#99FF99", alpha = 0.5) + 
          facet_wrap(~tw_window, nrow = 2) + 
          labs(title = "Summary of All Time Windows - Density Plots (min_est_over_se >= 5)") + 
          theme_minimal())
}

dev.off()
cat("Successfully generated PDF report:", output_pdf, "\n")
