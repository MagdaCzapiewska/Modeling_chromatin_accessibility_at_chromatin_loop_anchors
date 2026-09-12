#!/usr/bin/env Rscript

library(data.table)
library(yaml)

config_path <- file.path("config", "config.yml")
if (!file.exists(config_path)) {
  stop("File not found: config/config.yml")
}

config <- yaml::yaml.load_file(config_path)
datadir <- config$paths$datadir
resultsdir <- config$paths$resultsdir

loops_file <- file.path(datadir, "long_and_short_range_loops_D_mel.tsv")
counts_dir <- file.path(resultsdir, "counts", "counts_in_anchors")

if (!file.exists(loops_file)) {
  stop("File not found: ", loops_file)
}

if (!dir.exists(counts_dir)) {
  stop("Directory not found: ", counts_dir)
}

loops_table <- fread(loops_file)
if (!"loop_id" %in% names(loops_table)) {
  stop("Column 'loop_id' not found.")
}
loop_ids <- loops_table$loop_id

time_windows <- c("00-02", "02-04", "04-06", "06-08", "08-10", 
                  "10-12", "12-14", "14-16", "16-18", "18-20")

total_cells <- 0
sum_total_reads <- 0
sum_x_A1 <- 0
sum_x_A2 <- 0
zeros_x_A1 <- 0
zeros_x_A2 <- 0
files_processed <- 0


for (tw in time_windows) {
  for (loop_id in loop_ids) {
    

    fpath <- file.path(counts_dir, paste0("reads_", tw, "_", loop_id, ".tsv.gz"))
    if (!file.exists(fpath)) {
      fpath <- file.path(counts_dir, paste0("reads_", tw, "_", loop_id, ".tsv"))
    }
    
    if (!file.exists(fpath)) next
    
    dt <- fread(fpath, select = c("x_A1", "x_A2", "total_reads"))
    if (nrow(dt) == 0) next
    
    n <- nrow(dt)
    total_cells <- total_cells + n
    
    sum_total_reads <- sum_total_reads + sum(dt$total_reads, na.rm = TRUE)
    sum_x_A1        <- sum_x_A1 + sum(dt$x_A1, na.rm = TRUE)
    sum_x_A2        <- sum_x_A2 + sum(dt$x_A2, na.rm = TRUE)
    
    zeros_x_A1      <- zeros_x_A1 + sum(dt$x_A1 == 0, na.rm = TRUE)
    zeros_x_A2      <- zeros_x_A2 + sum(dt$x_A2 == 0, na.rm = TRUE)
    
    files_processed <- files_processed + 1
  }
}

if (total_cells == 0) {
  stop("No data found.")
}

mean_total_reads <- sum_total_reads / total_cells
mean_x_A1        <- sum_x_A1 / total_cells
mean_x_A2        <- sum_x_A2 / total_cells

pct_zeros_x_A1   <- (zeros_x_A1 / total_cells) * 100
pct_zeros_x_A2   <- (zeros_x_A2 / total_cells) * 100

# 6. Wypisanie wyników na konsolę
cat("\n==================================================\n")
cat("                   SUMMARY                          \n")
cat("==================================================\n")
cat(sprintf("Files:           %s\n", format(files_processed, big.mark = ",")))
cat(sprintf("Observations (N):    %s\n", format(total_cells, big.mark = ",")))
cat("--------------------------------------------------\n")
cat(sprintf("Mean total_reads:            %.2f\n", mean_total_reads))
cat(sprintf("Mean reads x_A1:   %.2f\n", mean_x_A1))
cat(sprintf("Mean reads x_A2:   %.2f\n", mean_x_A2))
cat("--------------------------------------------------\n")
cat(sprintf("Percentage of zeros in x_A1:              %.2f%%\n", pct_zeros_x_A1))
cat(sprintf("Percentage of zeros in x_A2:              %.2f%%\n", pct_zeros_x_A2))
cat("==================================================\n\n")
