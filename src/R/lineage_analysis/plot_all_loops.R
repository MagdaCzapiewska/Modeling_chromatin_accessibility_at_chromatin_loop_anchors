library(data.table)
library(ggplot2)
library(patchwork)
library(stringr)

config <- yaml::yaml.load_file(file.path("config", "config.yml"))
datadir <- config$paths$datadir
resultsdir <- config$paths$resultsdir

loops_file <- file.path(datadir, "long_and_short_range_loops_D_mel.tsv")
loops_data <- fread(loops_file)

root_to_leaf_paths <- file.path(resultsdir, "lineage_analysis", "root_to_leaf_paths.txt")
output_file <- file.path(resultsdir, "lineage_analysis", "correlation_by_lineage.pdf")

correlation_dir <- file.path(resultsdir, "MGLMfit_GDM_cor", "real_data", "pops", "init_1e-6")

time_windows <- c("00-02", "02-04", "04-06", "06-08", "08-10", "10-12", "12-14", "14-16", "16-18", "18-20")

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

all_cor <- all_cor[loop_id %in% loops_data$loop_id]
all_cor <- all_cor[min_est_over_se >= 5]

global_y_min <- min(all_cor$spearman_rho, na.rm = TRUE)
global_y_max <- max(all_cor$spearman_rho, na.rm = TRUE)

paths <- readLines(root_to_leaf_paths)
paths <- paths[paths != ""]

pdf(output_file, width = 14, height = 8)

for (p_idx in seq_along(paths)) {
  path_str <- paths[p_idx]
  
  path_nodes <- strsplit(path_str, "\\|")[[1]]
  path_nodes <- gsub("^\\(|\\)$", "", path_nodes)
  
  path_list <- lapply(path_nodes, function(node) {
    parts <- strsplit(node, ";")[[1]]
    data.table(time_window = parts[1], population = parts[2], lineage = parts[3])
  })
  path_dt <- rbindlist(path_list)
  
  cor_path <- merge(path_dt, all_cor, by = c("time_window", "population"), all.x = TRUE)
  counts <- cor_path[!is.na(spearman_rho), .(n_points = .N), by = time_window]
  

  template_dt <- data.table(time_window = time_windows)
  path_full <- merge(template_dt, path_dt, by = "time_window", all.x = TRUE)
  path_full <- merge(path_full, counts, by = "time_window", all.x = TRUE)
  path_full[is.na(n_points), n_points := 0]
  
  path_full[, x_label := ifelse(!is.na(population),
                                paste0(time_window, "\n", population, "\n", lineage, "\n(n=", n_points, ")"),
                                paste0(time_window, "\n-\n-\n(n=0)"))]
  
  path_full$x_label <- factor(path_full$x_label, levels = path_full$x_label)
  
  cor_plot <- merge(path_full[, .(time_window, x_label)], cor_path, by = "time_window", all.x = TRUE)
  
  has_data <- nrow(cor_plot[!is.na(spearman_rho)]) > 0
  
  p <- ggplot(cor_plot, aes(x = x_label, y = spearman_rho)) +
    scale_x_discrete(drop = FALSE) +
    coord_cartesian(ylim = c(global_y_min, global_y_max)) +
    labs(
      title = paste("Path", p_idx),
      x = "Time window \n Population \n Lineage \n Number of loops (n)",
      y = "Spearman Rho (min_est_over_se >= 5)"
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1, lineheight = 1.1),
      plot.subtitle = element_text(size = 8, color = "gray40")
    )
  
  if (has_data) {
    p <- p + 
      geom_violin(scale = "count", fill = "skyblue", color = "black", alpha = 0.6, na.rm = TRUE) +
      geom_jitter(position = position_jitter(width = 0.1, seed = 42), 
            alpha = 0.3, size = 0.8, color = "darkblue", na.rm = TRUE)
  } else {
    p <- p + 
      annotate("text", x = 5.5, y = mean(c(global_y_min, global_y_max)), 
               label = "No loops with min_est_over_se >= 5", color = "red")
  }
  
  print(p)
}

dev.off()
