library(data.table)
library(ggplot2)
library(patchwork)

if (!requireNamespace("yaml", quietly = TRUE)) install.packages("yaml") 

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript plot_naive_real.R <mode> <input_dir> <output_pdf>")
}

mode_arg   <- args[1]  # "all" or "pops"
input_dir  <- args[2]  
output_pdf <- args[3]  

tw_list <- c("00-02", "02-04", "04-06", "06-08", "08-10", "10-12", "12-14", "14-16", "16-18", "18-20")

config_path <- file.path("config", "config.yml")
if (file.exists(config_path)) {
  config <- yaml::yaml.load_file(config_path)
  resultsdir <- config$paths$resultsdir
  cardinality_file <- file.path(resultsdir, "cluster_cardinality.tsv.gz")
} else {
  cardinality_file <- "results/cluster_cardinality.tsv.gz"
}

plot_naive_pair <- function(dt, title_prefix = "", text_size_factor = 1, is_pop_mode = FALSE, total_cells = NULL) {
  n_loops <- nrow(dt[!is.na(spearman_rho)])
  cell_info <- if(!is.null(total_cells)) paste0(", cells = ", total_cells) else ""
  
  spearman_title <- if(is_pop_mode) paste0(title_prefix, "\n[Spearman] (n=", n_loops, cell_info, ")") else "Spearman Rho"
  pearson_title  <- if(is_pop_mode) paste0(title_prefix, "\n[Pearson]") else "Pearson R"
  
  p_spearman <- ggplot(dt, aes(x = spearman_rho)) +
    geom_histogram(binwidth = 0.05, boundary = 0, fill = "#708090", color = "white", linewidth = 0.1) +
    xlim(-1, 1) +
    labs(title = spearman_title, x = "Rho", y = "Count") +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 8 * text_size_factor, face = if(is_pop_mode) "bold" else "plain"), 
      axis.title = element_text(size = 7 * text_size_factor)
    )
  
  p_pearson <- ggplot(dt, aes(x = pearson_r)) +
    geom_histogram(binwidth = 0.05, boundary = 0, fill = "#4682B4", color = "white", linewidth = 0.1) +
    xlim(-1, 1) +
    labs(title = pearson_title, x = "R", y = "Count") +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 8 * text_size_factor, face = if(is_pop_mode) "bold" else "plain"), 
      axis.title = element_text(size = 7 * text_size_factor), 
      axis.title.y = element_blank()
    )
  
  if (is_pop_mode) {
    return(p_spearman + p_pearson)
  } else {
    full_title <- paste0(title_prefix, " (loops = ", n_loops, cell_info, ")")
    return(p_spearman + p_pearson + plot_annotation(title = full_title, theme = theme(plot.title = element_text(size = 12, face = "bold"))))
  }
}

if(!dir.exists(dirname(output_pdf))) dir.create(dirname(output_pdf), recursive = TRUE)

has_cardinality <- FALSE
if (file.exists(cardinality_file)) {
  dt_cardinality <- fread(cardinality_file)
  cluster_col <- grep("cluster", names(dt_cardinality), value = TRUE, ignore.case = TRUE)[1]
  
  if (length(cluster_col) > 0) {
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
  warning(paste("Warning: no cardinality file:", cardinality_file))
}

all_data_list <- list()
pdf(output_pdf, width = 16, height = 12) 

for (tw in tw_list) {
  path <- file.path(input_dir, paste0("cor_", tw, ".tsv.gz"))
  if(!file.exists(path)) next
  
  dt <- fread(path)
  dt[, tw_window := tw]
  all_data_list[[tw]] <- dt
  
  global_cells <- NULL
  if (mode_arg == "all" && has_cardinality) {
    global_cells <- sum(dt_card_long[tw_window == tw, cells], na.rm = TRUE)
  }
  
  if (mode_arg == "all") {
    print(plot_naive_pair(dt, paste("Naive Real TW:", tw), text_size_factor = 1.5, is_pop_mode = FALSE, total_cells = global_cells))
  } else {
    if (!"population" %in% names(dt) || nrow(dt) == 0) next
    dt[, pop_id := as.integer(sub("^([0-9]+)_.*", "\\1", population))]
    
    pop_info <- unique(dt[!is.na(population), .(pop_id, population)])
    
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
    
    p_list <- list()
    for(i in seq_len(nrow(pop_info))) {
      pid   <- pop_info$pop_id[i]
      pname <- pop_info$population[i] 
      pcell <- if(has_cardinality) pop_info$cells[i] else NULL
      sub_dt <- dt[pop_id == pid]
      
      p_list[[i]] <- plot_naive_pair(sub_dt, title_prefix = pname, text_size_factor = 1.1, is_pop_mode = TRUE, total_cells = pcell)
    }
    
    pops_per_page <- 8  
    plot_chunks <- split(p_list, ceiling(seq_along(p_list) / pops_per_page))
    
    for (page_idx in seq_along(plot_chunks)) {
      page_plots <- plot_chunks[[page_idx]]
      
      print(wrap_plots(page_plots, ncol = 2) + 
            plot_layout(guides = "collect") +
            plot_annotation(
              title = paste0("Time Window: ", tw, " - Naive Correlations per Population"),
              subtitle = paste0("Page ", page_idx, " of ", length(plot_chunks), " for this Time Window"),
              theme = theme(
                plot.title = element_text(size = 15, face = "bold"),
                plot.subtitle = element_text(size = 11, face = "italic")
              )
            ))
    }
  }
}

if (mode_arg == "all" && length(all_data_list) > 0) {
  big_dt <- rbindlist(all_data_list)
  
  p_all_s <- ggplot(big_dt, aes(x = spearman_rho)) +
    geom_histogram(binwidth = 0.05, boundary = 0, fill = "#708090", color = "white", linewidth = 0.1) +
    facet_wrap(~tw_window, nrow = 2) +
    labs(title = "Summary: Naive Spearman Rho", x = "Rho", y = "Count") + 
    theme_minimal()
  
  p_all_p <- ggplot(big_dt, aes(x = pearson_r)) +
    geom_histogram(binwidth = 0.05, boundary = 0, fill = "#4682B4", color = "white", linewidth = 0.1) +
    facet_wrap(~tw_window, nrow = 2) +
    labs(title = "Summary: Naive Pearson R", x = "R", y = "Count") + 
    theme_minimal()
  
  print(p_all_s)
  print(p_all_p)
}

dev.off()
cat("Successfully generated Naive Real PDF report with split population pages and cell tracking:", output_pdf, "\n")
