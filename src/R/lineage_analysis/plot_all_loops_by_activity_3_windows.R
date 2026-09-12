library(data.table)
library(ggplot2)
library(stringr)

config <- yaml::yaml.load_file(file.path("config", "config.yml"))
datadir <- config$paths$datadir
resultsdir <- config$paths$resultsdir

loops_file <- file.path(datadir, "long_and_short_range_loops_D_mel.tsv")
loops_data <- fread(loops_file)

root_to_leaf_paths <- file.path(resultsdir, "lineage_analysis", "root_to_leaf_paths.txt")
output_file <- file.path(resultsdir, "lineage_analysis", "correlation_by_lineage_with_union_activity_3_windows.pdf")
correlation_dir <- file.path(resultsdir, "MGLMfit_GDM_cor", "real_data", "pops", "init_1e-6")

time_windows <- c("06-08", "10-12", "14-16")

act_06_08 <- loops_data[, .(
  loop_id, 
  time_window = "06-08", 
  is_active = as.integer((`Dmel_6-8h_Neuroblasts` == 1) | (`Dmel_6-8h_Neurons` == 1) | (`Dmel_6-8h_Glia` == 1))
)]

act_10_12 <- loops_data[, .(
  loop_id, 
  time_window = "10-12", 
  is_active = as.integer((`Dmel_10-12h_Neuroblasts` == 1) | (`Dmel_10-12h_Neurons` == 1) | (`Dmel_10-12h_Glia` == 1))
)]

act_14_16 <- loops_data[, .(
  loop_id, 
  time_window = "14-16", 
  is_active = as.integer((`Dmel_14-16h_Neuroblasts` == 1) | (`Dmel_14-16h_Neurons` == 1) | (`Dmel_14-16h_Glia` == 1))
)]

activity_map <- rbindlist(list(act_06_08, act_10_12, act_14_16))

cor_list <- list()
for (tw in time_windows) {
  file_path <- file.path(correlation_dir, paste0("cor_", tw, ".tsv.gz"))
  if (file.exists(file_path)) {
    dt <- fread(file_path)
    cor_list[[tw]] <- dt
  } else {
    warning(paste("File not exists:", file_path))
  }
}
all_cor <- rbindlist(cor_list)

all_cor <- all_cor[loop_id %in% loops_data$loop_id]
all_cor <- all_cor[min_est_over_se >= 5]
all_cor <- merge(all_cor, activity_map, by = c("loop_id", "time_window"), all.x = TRUE)

global_y_min <- min(all_cor$spearman_rho, na.rm = TRUE)
global_y_max <- max(all_cor$spearman_rho, na.rm = TRUE)

paths <- readLines(root_to_leaf_paths)
paths <- paths[paths != ""]

template_dt <- data.table(time_window = time_windows, x_base = 1:3)

pdf(output_file, width = 10, height = 7)

for (p_idx in seq_along(paths)) {
  path_str <- paths[p_idx]
  
  path_nodes <- strsplit(path_str, "\\|")[[1]]
  path_nodes <- gsub("^\\(|\\)$", "", path_nodes)
  
  path_list <- lapply(path_nodes, function(node) {
    parts <- strsplit(node, ";")[[1]]
    data.table(time_window = parts[1], population = parts[2], lineage = parts[3])
  })
  path_dt <- rbindlist(path_list)
  
  path_dt <- path_dt[time_window %in% time_windows]
  
  cor_path <- merge(path_dt, all_cor, by = c("time_window", "population"), all.x = TRUE)
  
  plot_data <- copy(cor_path)
  plot_data[, split_var := factor(is_active, levels = c("0", "1"))]
  

  counts_split <- plot_data[!is.na(spearman_rho), .(n = .N), by = .(time_window, split_var)]
  counts_wide <- dcast(counts_split, time_window ~ split_var, value.var = "n", fill = 0, drop = FALSE)
  
  if (!"0" %in% names(counts_wide)) counts_wide[, `0` := 0]
  if (!"1" %in% names(counts_wide)) counts_wide[, `1` := 0]
  setnames(counts_wide, c("0", "1"), c("n_0", "n_1"))
  
  path_full <- merge(template_dt, path_dt, by = "time_window", all.x = TRUE)
  path_full <- merge(path_full, counts_wide, by = "time_window", all.x = TRUE)
  path_full[is.na(n_0), n_0 := 0]
  path_full[is.na(n_1), n_1 := 0]
  
  path_full[, x_label := ifelse(!is.na(population),
                                paste0(time_window, "\n", population, "\n", lineage, "\n(0: ", n_0, " | 1: ", n_1, ")"),
                                paste0(time_window, "\n-\n-\n(0: 0 | 1: 0)"))]
  
  path_full <- path_full[order(x_base)]
  path_full$x_label <- factor(path_full$x_label, levels = path_full$x_label)
  
  plot_data <- merge(plot_data, path_full[, .(time_window, x_label)], by = "time_window", all.x = TRUE)
  
  has_data <- nrow(plot_data[!is.na(spearman_rho) & !is.na(split_var)]) > 0
  
  p <- ggplot(plot_data, aes(x = x_label, y = spearman_rho, fill = split_var, color = split_var)) +
    scale_x_discrete(limits = path_full$x_label, drop = FALSE) +
    coord_cartesian(ylim = c(global_y_min, global_y_max)) +
    scale_fill_manual(values = c("0" = "skyblue", "1" = "salmon"), drop = FALSE) +
    scale_color_manual(values = c("0" = "darkblue", "1" = "darkred"), drop = FALSE) +
    labs(
      title = paste("Path", p_idx),
      subtitle = "Loop is 'Active' (1) if it appears in Neuroblasts, Neurons, OR Glia for that specific time window.",
      x = "Time window \n Population \n Lineage \n (0: Inactive loops | 1: Active loops)",
      y = "Spearman Rho",
      fill = "Union Activity\n(Any tissue)",
      color = "Union Activity\n(Any tissue)"
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1, lineheight = 1.1),
      legend.position = "right"
    )
  
  if (has_data) {
    p <- p + 
      geom_violin(position = position_dodge(width = 0.8), scale = "count", 
                  alpha = 0.6, color = "black", na.rm = TRUE) +
      geom_jitter(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8, seed = 42), 
                  alpha = 0.4, size = 0.8, na.rm = TRUE)
  } else {
    p <- p + annotate("text", x = 2, y = mean(c(global_y_min, global_y_max)), 
                      label = "No valid loops in these windows", color = "red")
  }
  
  print(p)
}

dev.off()
