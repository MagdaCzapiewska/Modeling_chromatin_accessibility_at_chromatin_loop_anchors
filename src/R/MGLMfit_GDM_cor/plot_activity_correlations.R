library(data.table)
library(ggplot2)
library(patchwork)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage: Rscript plot_activity_correlations.R <mode> <init> <input_dir> <output_pdf>")
}

mode_arg   <- args[1]  # "all" or "pops"
init_val   <- args[2]  # "default" or "1e-6"
input_dir  <- args[3]  # directory with files cor_{tw}.tsv.gz
output_pdf <- args[4]  # output PDF path

config <- yaml::yaml.load_file(file.path("config", "config.yml"))
datadir <- config$paths$datadir

tw_list <- c("06-08", "10-12", "14-16")

loops_file <- file.path(datadir, "long_and_short_range_loops_D_mel.tsv")
loops_data <- fread(loops_file)

color_map <- c(
  "10" = "#99FF99",
  "5"  = "#FFFF99",
  "2"  = "#FFCC99",
  "1"  = "#FFB3B3",
  "0"  = "#99CCFF",
  "NA" = "#D3D3D3"
)

act_color   <- "#99FF99"
inact_color <- "#FFCC99"

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

plot_type_1 <- function(dt, title_prefix = "") {
  n_loops <- nrow(dt[!is.na(spearman_rho)])
  full_title <- paste0(title_prefix, "\n(loops = ", n_loops, ")")
  
  ggplot(dt, aes(x = spearman_rho, fill = color_group)) +
    geom_histogram(binwidth = 0.05, boundary = 0, color = "white", linewidth = 0.1) +
    scale_fill_manual(values = color_map, drop = FALSE, name = "min_est_over_se") +
    xlim(-1, 1) +
    labs(title = full_title, x = "Spearman Rho", y = "Count") +
    theme_minimal() +
    theme(plot.title = element_text(size = 11), axis.title = element_text(size = 7))
}

plot_type_2 <- function(dt, title_prefix = "", fill_color = "#99FF99") {
  dt_sub <- dt[as.numeric(as.character(color_group)) >= 5 & !is.na(as.numeric(as.character(color_group)))]
  n_loops <- nrow(dt_sub[!is.na(spearman_rho)])
  full_title <- paste0(title_prefix, "\n(loops >= 5: ", n_loops, ")")
  
  p <- ggplot(dt_sub, aes(x = spearman_rho)) +
    geom_density(fill = fill_color, alpha = 0.5) +
    xlim(-1, 1) +
    labs(title = full_title, subtitle = "min_est_over_se >= 5", x = "Spearman Rho", y = "Density") +
    theme_minimal() +
    theme(plot.title = element_text(size = 11), axis.title = element_text(size = 7))
  
  if(nrow(dt_sub) < 2) {
    p <- p + annotate("text", x = 0, y = 0.5, label = "No data (>=5)", size = 2)
  }
  return(p)
}


if(!dir.exists(dirname(output_pdf))) dir.create(dirname(output_pdf), recursive = TRUE)

categories <- list(
  "Anywhere"    = "loop_active_anywhere",
  "Neurons"     = "loop_active_neur",
  "Neuroblasts" = "loop_active_nb",
  "Glia"        = "loop_active_glia"
)

pdf(output_pdf, width = 18, height = 12)

for (tw in tw_list) {
  path <- file.path(input_dir, paste0("cor_", tw, ".tsv.gz"))
  if(!file.exists(path)) next
  
  dt <- fread(path)
  dt <- process_data(dt)
  
  tw_label <- paste0(as.numeric(substr(tw, 1, 2)), "-", as.numeric(substr(tw, 4, 5)), "h")
  col_nb   <- paste0("Dmel_", tw_label, "_Neuroblasts")
  col_neur <- paste0("Dmel_", tw_label, "_Neurons")
  col_glia <- paste0("Dmel_", tw_label, "_Glia")
  
  if (!all(c(col_nb, col_neur, col_glia) %in% names(loops_data))) next
  
  loops_sub <- loops_data[, .(loop_id, 
                              loop_active_nb = get(col_nb), 
                              loop_active_neur = get(col_neur), 
                              loop_active_glia = get(col_glia))]
  loops_sub[, loop_active_anywhere := as.integer(loop_active_nb | loop_active_neur | loop_active_glia)]
  
  dt <- merge(dt, loops_sub, by = "loop_id", all.x = TRUE)
  
  for (cat_name in names(categories)) {
    cat_col <- categories[[cat_name]]
    dt_active <- dt[get(cat_col) == 1]
    dt_inactive <- dt[get(cat_col) == 0]
    
    if (mode_arg == "all") {
      p1_act <- plot_type_1(dt_active, paste(cat_name, "Active"))
      p1_inact <- plot_type_1(dt_inactive, paste(cat_name, "Inactive"))
      print(p1_act + p1_inact + plot_annotation(title = paste("TW:", tw, "| Cat:", cat_name, "- Histograms")))
      
      p2_act <- plot_type_2(dt_active, paste(cat_name, "Active"), fill_color = act_color)
      p2_inact <- plot_type_2(dt_inactive, paste(cat_name, "Inactive"), fill_color = inact_color)
      
      dt_comp <- dt[get(cat_col) %in% c(0, 1) & as.numeric(as.character(color_group)) >= 5]
      p2_overlay <- ggplot(dt_comp, aes(x = spearman_rho, fill = factor(get(cat_col)))) +
        geom_density(alpha = 0.4) + xlim(-1, 1) + theme_minimal() +
        scale_fill_manual(values = c("0" = inact_color, "1" = act_color), labels = c("Inactive", "Active"), name = "Status") +
        labs(title = "Overlay Comparison", subtitle = "min_est_over_se >= 5")
      
      print(p2_act + p2_inact + p2_overlay + plot_layout(ncol = 3) + 
              plot_annotation(title = paste("TW:", tw, "| Cat:", cat_name, "- Density Plots")))
      
    } else {
      pop_info <- unique(dt[!is.na(population), .(pop_id, pop_name)])
      setorder(pop_info, pop_id)
      
      if(nrow(pop_info) == 0) next
      
      hist_list <- list()
      dens_list <- list()
      
      for (i in seq_len(nrow(pop_info))) {
        pid <- pop_info$pop_id[i]; pname <- pop_info$pop_name[i]
        sub_act <- dt[pop_id == pid & get(cat_col) == 1]
        sub_inact <- dt[pop_id == pid & get(cat_col) == 0]
        
        hist_list[[length(hist_list) + 1]] <- plot_type_1(sub_act, paste(pname, "(Act)")) + theme(legend.position = "none")
        hist_list[[length(hist_list) + 1]] <- plot_type_1(sub_inact, paste(pname, "(Inact)")) + theme(legend.position = "none")
        
        dens_list[[length(dens_list) + 1]] <- plot_type_2(sub_act, paste(pname, "(Act)"), fill_color = act_color)
        dens_list[[length(dens_list) + 1]] <- plot_type_2(sub_inact, paste(pname, "(Inact)"), fill_color = inact_color)
        
        sub_comp <- dt[pop_id == pid & get(cat_col) %in% c(0, 1) & as.numeric(as.character(color_group)) >= 5]
        p_pop_overlay <- ggplot(sub_comp, aes(x = spearman_rho, fill = factor(get(cat_col)))) +
          geom_density(alpha = 0.4) + xlim(-1, 1) + theme_minimal() +
          scale_fill_manual(values = c("0" = inact_color, "1" = act_color), guide = "none") +
          labs(title = paste(pname, "(Comp)"), x = "Rho") + theme(plot.title = element_text(size = 11))
        
        dens_list[[length(dens_list) + 1]] <- p_pop_overlay
      }
      
      hist_pages <- split(hist_list, ceiling(seq_along(hist_list) / 8))
      for(page in hist_pages) {
        print(wrap_plots(page, ncol = 4, nrow = 2) + 
                plot_annotation(title = paste("TW:", tw, "| Cat:", cat_name, "- Histograms")))
      }
      
      dens_pages <- split(dens_list, ceiling(seq_along(dens_list) / 6))
      for(page in dens_pages) {
        print(wrap_plots(page, ncol = 3, nrow = 2) + 
                plot_annotation(title = paste("TW:", tw, "| Cat:", cat_name, "- Density & Overlays")))
      }
    }
  }
}
dev.off()
cat("Successfully generated activity report:", output_pdf, "\n")
