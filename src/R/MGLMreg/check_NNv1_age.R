library(data.table)

config_path <- file.path("config", "config.yml")
config <- yaml::yaml.load_file(config_path)

datadir <- config$paths$datadir
resultsdir <- config$paths$resultsdir

prev_dir <- file.path("..", "Genomics", "results", "dmel", "embrio", "GDM_reads_with_NNv1_age")
current_dir <- file.path(resultsdir, "counts", "counts_in_anchors")

loops_file <- file.path(datadir, "long_and_short_range_loops_D_mel.tsv")
loops_table <- fread(loops_file)
loops_of_interest <- loops_table$loop_id

time_windows <- c("00-02", "02-04", "04-06", "06-08", "08-10", "10-12", "12-14", "14-16", "16-18", "18-20")

total_checked <- 0
mismatches <- 0
missing_files <- 0

cat("=== Starting validation ===\n")

for (tw in time_windows) {
  for (loop_id in loops_of_interest) {
    
    file_prev <- file.path(prev_dir, paste0("reads_age_", tw, "_", loop_id, ".tsv"))
    file_curr <- file.path(current_dir, paste0("reads_", tw, "_", loop_id, ".tsv.gz"))
    
    if (!file.exists(file_prev) || !file.exists(file_curr)) {
      if (!file.exists(file_prev) && !file.exists(file_curr)) next

      cat(sprintf("No pair of files for TW: %s, Loop: %s (Prev exists: %s, Curr exists: %s)\n", 
                  tw, loop_id, file.exists(file_prev), file.exists(file_curr)))
      missing_files <- missing_files + 1
      next
    }
    
    total_checked <- total_checked + 1
    
    dt_prev <- fread(file_prev)
    dt_curr <- fread(file_curr)

    if ("cell" %in% names(dt_prev)) {
      setnames(dt_prev, "cell", "barcode")
    }
    
    if (!all(names(dt_prev) %in% names(dt_curr)) || !all(names(dt_curr) %in% names(dt_prev))) {
      cat(sprintf("Unrecoverable difference in column names for TW: %s, Loop: %s\n", tw, loop_id))
      cat("   Prev (standardized):", paste(names(dt_prev), collapse=", "), "\n")
      cat("   Curr:", paste(names(dt_curr), collapse=", "), "\n")
      mismatches <- mismatches + 1
      next
    }
    
    sorted_cols <- sort(names(dt_curr))
    setcolorder(dt_prev, sorted_cols)
    setcolorder(dt_curr, sorted_cols)
    
    diff_check <- all.equal(dt_prev, dt_curr, check.attributes = FALSE)
    
    if (!isTRUE(diff_check)) {
      cat(sprintf("Difference in content for TW: %s, Loop: %s\n", tw, loop_id))
      cat("   Details:", head(diff_check, 1), "\n")
      mismatches <- mismatches + 1
    }
  }
}

cat("\n=== Summary ===\n")
cat("Number of pairs: ", total_checked, "\n")
cat("Files with mismatches:  ", mismatches, "\n")
cat("Not paired files:   ", missing_files, "\n")

if (mismatches == 0 && missing_files == 0 && total_checked > 0) {
  cat("Success! All files are aligned and contain identical data.\n")
} else {
  cat("Failure.\n")
}
