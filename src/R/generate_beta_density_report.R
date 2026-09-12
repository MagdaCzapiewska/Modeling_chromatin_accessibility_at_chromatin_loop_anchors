library(data.table)
library(dplyr)
library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop("Missing required arguments!\nUsage: Rscript generate_beta_density_report.R <loop_id> <xmin> <xmax>")
}

loop_id <- args[1]
xmin    <- as.numeric(args[2])
xmax    <- as.numeric(args[3])

if (is.na(xmin) || is.na(xmax)) {
  stop("Error: <xmin> and <xmax> must be numeric values!")
}
if (xmin >= xmax) {
  stop("Error: <xmin> must be strictly less than <xmax>!")
}

config_path <- file.path("config", "config.yml")
config <- yaml::yaml.load_file(config_path)

resultsdir <- config$paths$resultsdir
beta_dir   <- file.path(resultsdir, "parameters", "beta_distribution_of_counts_in_anchors")
beta_plots_dir <- file.path(resultsdir, "parameters", "beta_distribution_of_counts_in_anchors_plots")
if (!dir.exists(beta_plots_dir)) {
  dir.create(beta_plots_dir, recursive = TRUE, showWarnings = FALSE)
}
output_pdf <- file.path(beta_plots_dir, paste0("report_beta_densities_loop_", loop_id, ".pdf"))

time_windows <- c("00-02", "02-04", "04-06", "06-08", "08-10", 
                  "10-12", "12-14", "14-16", "16-18", "18-20")

x_seq <- seq(xmin, xmax, length.out = 500)

message("=== Starting Beta Distribution PDF Report Generation ===")
message("Target Loop ID: ", loop_id)
message("Input directory: ", beta_dir)

pdf(output_pdf, width = 14, height = 8.5)

# -------------------------------------------------------------------------
# PAGES 1 & 2: GLOBAL OVERVIEWS ACROSS ALL TIME WINDOWS (POPULATION 'ALL')
# -------------------------------------------------------------------------
for (target_anchor in c("A1", "A2")) {
  message("Generating global overview page for Anchor: ", target_anchor, "...")
  plot_data_p12 <- list()
  
  for (tw in time_windows) {
    file_path <- file.path(beta_dir, paste0("beta_parameters_", tw, "_", loop_id, "_all_and_by_population.tsv.gz"))
    if (!file.exists(file_path)) next
    
    dt <- fread(file_path)[population == "all" & anchor == target_anchor]
    if (nrow(dt) == 0 || is.na(dt$alpha_prior) || is.na(dt$beta_prior)) next
    
    current_alpha <- dt$alpha_prior
    current_beta  <- dt$beta_prior
    
    mu <- current_alpha / (current_alpha + current_beta)
    nu <- current_alpha + current_beta
    
    y_vals <- dbeta(x_seq, shape1 = current_alpha, shape2 = current_beta)
    
    lbl <- paste0(tw, " (a=", round(current_alpha, 2), ", b=", round(current_beta, 2), "\nm=", round(mu, 6), ", n=", round(nu, 1), ")\n")
    plot_data_p12[[tw]] <- data.frame(x = x_seq, y = y_vals, tw = tw, legend_lbl = lbl)
  }
  
  if (length(plot_data_p12) > 0) {
    df_p12 <- bind_rows(plot_data_p12)
    
    p_global <- ggplot(df_p12, aes(x = x, y = y, color = legend_lbl)) +
      geom_line(linewidth = 1.2, alpha = 0.85) +
      labs(
        title = paste("Beta Prior Distributions Across All Time Windows - Anchor", target_anchor),
        subtitle = paste("Global population ('all') density curves for Loop:", loop_id),
        x = "Probability of Success (p)",
        y = "Probability Density",
        color = "Time Window (Parameters)"
      ) +
      theme_minimal(base_size = 12) +
      theme(
        plot.title = element_text(face = "bold", size = 16, color = "#0f172a"),
        plot.subtitle = element_text(size = 11, color = "#475569"),
        legend.position = "right",
        legend.text = element_text(size = 9.5),
        panel.grid.minor = element_blank()
      ) +
      scale_color_brewer(palette = "Spectral")
    
    print(p_global)
  } else {
    message("WARNING: No data available for Anchor ", target_anchor, " global overview page.")
  }
}

# -------------------------------------------------------------------------
# PAGES 3+: PER-TIME-WINDOW DECOMPOSITION (2 PAGES PER WINDOW: A1 THEN A2)
# -------------------------------------------------------------------------
for (tw in time_windows) {
  file_path <- file.path(beta_dir, paste0("beta_parameters_", tw, "_", loop_id, "_all_and_by_population.tsv.gz"))
  if (!file.exists(file_path)) {
    message("Skipped window ", tw, " - .tsv.gz parameter file not found")
    next
  }
  
  dt_tw_all <- fread(file_path)
  message("Processing window ", tw, "...")

  for (target_anchor in c("A1", "A2")) {
    dt_anchor <- dt_tw_all[anchor == target_anchor]
    if (nrow(dt_anchor) == 0) next
    
    plot_data_tw <- list()
    for (i in seq_len(nrow(dt_anchor))) {
      row <- dt_anchor[i]
      if (is.na(row$alpha_prior) || is.na(row$beta_prior)) next
      
      sub_alpha <- row$alpha_prior
      sub_beta  <- row$beta_prior
      
      mu_pop <- sub_alpha / (sub_alpha + sub_beta)
      nu_pop <- sub_alpha + sub_beta
      
      y_vals_sub <- dbeta(x_seq, shape1 = sub_alpha, shape2 = sub_beta)
      
      lbl <- paste0(row$population, " (a=", round(sub_alpha, 2), ", b=", round(sub_beta, 2), "\nm=", round(mu_pop, 6), ", n=", round(nu_pop, 1), ")\n")
      plot_data_tw[[row$population]] <- data.frame(x = x_seq, y = y_vals_sub, pop = row$population, legend_lbl = lbl)
    }
    
    if (length(plot_data_tw) == 0) next
    
    df_tw <- bind_rows(plot_data_tw)
    df_tw$line_weight <- ifelse(df_tw$pop == "all", "all_group", "subgroup")
    
    unique_lbls <- unique(df_tw$legend_lbl)
    
    num_subgroups <- length(unique_lbls) - 1
    subgroup_colors <- scales::hue_pal()(num_subgroups)
    
    custom_colors <- c()
    subgroup_idx <- 1
    
    for (lbl in unique_lbls) {
      if (grepl("^all \\(", lbl)) {
        custom_colors[lbl] <- "#000000"
      } else {
        custom_colors[lbl] <- subgroup_colors[subgroup_idx]
        subgroup_idx <- subgroup_idx + 1
      }
    }
    # -----------------------------------------------------------
    
    p_tw <- ggplot(df_tw, aes(x = x, y = y, color = legend_lbl, size = line_weight)) +
      geom_line(alpha = 0.85) +
      scale_size_manual(values = c("all_group" = 1.7, "subgroup" = 0.9), guide = "none") +
      scale_color_manual(values = custom_colors) + 
      labs(
        title = paste("Time Window:", tw, "— Anchor:", target_anchor),
        subtitle = paste("Comparison of fitted Beta prior distributions for Loop:", loop_id),
        x = "Probability of Success (p)",
        y = "Probability Density",
        color = "Population (Parameters)"
      ) +
      theme_minimal(base_size = 12) +
      theme(
        plot.title = element_text(face = "bold", size = 16, color = "#0f172a"),
        plot.subtitle = element_text(size = 11, color = "#475569"),
        legend.position = "right",
        legend.text = element_text(size = 9),
        panel.grid.minor = element_blank()
      )
    
    print(p_tw)
  }
}

dev.off()
message("=== Success: PDF Report generated successfully! ===")
message("Output file path: ", output_pdf)
