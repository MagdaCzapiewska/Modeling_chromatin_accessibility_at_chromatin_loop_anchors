library(data.table)
library(ggplot2)
library(stringr)

config <- yaml::yaml.load_file(file.path("config", "config.yml"))
datadir <- config$paths$datadir
resultsdir <- config$paths$resultsdir

loops_file <- file.path(datadir, "long_and_short_range_loops_D_mel.tsv")
parents_file <- file.path(resultsdir, "lineage_analysis", "lineage_parents_mapping.tsv")
output_file <- file.path(resultsdir, "lineage_analysis", "correlation_by_loop.pdf")
correlation_dir <- file.path(resultsdir, "MGLMfit_GDM_cor", "real_data", "pops", "init_1e-6")

time_windows <- c("00-02", "02-04", "04-06", "06-08", "08-10", "10-12", "12-14", "14-16", "16-18", "18-20")

loops_data <- fread(loops_file)

parents_mapping <- fread(parents_file)
parents_mapping[, node_id := paste(time_window, population, sep = "__")]
parents_mapping[, parent_node_id := ifelse(parent_time_window == "ROOT", 
                                           NA_character_, 
                                           paste(parent_time_window, parent_population, sep = "__"))]

parent_dict <- setNames(parents_mapping$parent_node_id, parents_mapping$node_id)

cor_list <- list()
for (tw in time_windows) {
  file_path <- file.path(correlation_dir, paste0("cor_", tw, ".tsv.gz"))
  if (file.exists(file_path)) {
    cor_list[[tw]] <- fread(file_path)
  } else {
    warning(paste("File not existing:", file_path))
  }
}
all_cor <- rbindlist(cor_list)

all_cor <- all_cor[loop_id %in% loops_data$loop_id]
all_cor <- all_cor[min_est_over_se >= 5]

all_cor[, node_id := paste(time_window, population, sep = "__")]
all_cor[, tissue := sub("^[^_]+_", "", population)]
all_cor[, x_pos := match(time_window, time_windows)]

set.seed(42)

all_cor[, x_jittered := x_pos + runif(.N, min = -0.15, max = 0.15)]
all_cor[, y_jittered := spearman_rho + runif(.N, min = -0.005, max = 0.005)]

global_y_min <- min(all_cor$y_jittered, na.rm = TRUE)
global_y_max <- max(all_cor$y_jittered, na.rm = TRUE)
y_margin <- (global_y_max - global_y_min) * 0.05
y_limits <- c(global_y_min - y_margin, global_y_max + y_margin)

loops_to_plot <- unique(all_cor$loop_id)
loops_to_plot <- loops_to_plot[order(as.numeric(gsub("[^0-9]", "", loops_to_plot)))]

pdf(output_file, width = 12, height = 8)

for (l_id in loops_to_plot) {

  message(l_id, "\n")
  
  loop_data <- all_cor[loop_id == l_id]
  valid_nodes <- loop_data$node_id

  edges <- data.table(child_node = character(), parent_node = character())
  
  for (n in valid_nodes) {
    curr_parent <- parent_dict[[n]]
    
    while (!is.na(curr_parent) && !(curr_parent %in% valid_nodes)) {
      curr_parent <- parent_dict[[curr_parent]]
    }
    
    if (!is.na(curr_parent)) {
      edges <- rbind(edges, data.table(child_node = n, parent_node = curr_parent))
    }
  }

  p <- ggplot() +
    scale_x_continuous(breaks = 1:length(time_windows), 
                       labels = time_windows, 
                       limits = c(0.8, length(time_windows) + 0.2)) +
    coord_cartesian(ylim = y_limits) +
    labs(
      title = paste("Correlation (Spearman) for loop:", l_id),
      subtitle = "min_est_over_se >= 5",
      x = "Time window",
      y = "Spearman Rho",
      color = "Tissue"
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "right"
    )
  
  if (nrow(edges) > 0) {

    edges_coords <- merge(edges, loop_data[, .(node_id, x_jittered, y_jittered)], 
                          by.x = "child_node", by.y = "node_id")
    setnames(edges_coords, c("x_jittered", "y_jittered"), c("xend", "yend"))
    
    edges_coords <- merge(edges_coords, loop_data[, .(node_id, x_jittered, y_jittered)], 
                          by.x = "parent_node", by.y = "node_id")
    setnames(edges_coords, c("x_jittered", "y_jittered"), c("x", "y"))
    
    p <- p + geom_segment(data = edges_coords, 
                          aes(x = x, y = y, xend = xend, yend = yend), 
                          color = "gray60", alpha = 0.6, linewidth = 0.5)
  }
  

  p <- p + geom_point(data = loop_data, 
                      aes(x = x_jittered, y = y_jittered, color = tissue), 
                      size = 3, alpha = 0.9)
  
  print(p)
}

dev.off()
