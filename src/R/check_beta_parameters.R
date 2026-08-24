library(data.table)
library(dplyr)

config_path <- file.path("config", "config.yml")
config <- yaml::yaml.load_file(config_path)

resultsdir <- config$paths$resultsdir

loops_ids <- paste0("L", 1:417)

new_results <- file.path(resultsdir, "parameters", "beta_distribution_of_counts_in_anchors")
old_results <- file.path("..", "Genomics", "results", "dmel", "embrio", "p_fitting_random_loops", "all_loops_joint_all_pop_neural", "processed_data")

time_windows <- c("00-02", "02-04", "04-06", "06-08", "08-10", "10-12", "12-14", "14-16", "16-18", "18-20")

expected_cols <- c("population", "anchor", "alpha_prior", "beta_prior", "n_cells")

for (tw in time_windows) {
  cat("\n--- Checking time window:", tw, "---\n")
  
  for (loop_id in loops_ids) {
    cat("Loop id: ", loop_id, "\n")
    
    old_file <- file.path(old_results, paste0("ebbinom_prior_", tw, "_", loop_id, "_by_population.tsv"))
    new_file <- file.path(new_results, paste0("beta_parameters_", tw, "_", loop_id, "_all_and_by_population.tsv.gz"))

    if (!file.exists(old_file)) {
      stop(paste("ERROR: File does not exist in old directory:", old_file))
    }
    if (!file.exists(new_file)) {
      stop(paste("ERROR: File does not exist in new directory:", new_file))
    }

    dt_old <- fread(old_file)
    dt_new <- fread(new_file)

    if (!all(names(dt_old) == expected_cols)) {
      stop(paste("ERROR: Invalid columns in old file:", old_file, 
                 "\nExpected:", paste(expected_cols, collapse=", "), 
                 "\nFound:", paste(names(dt_old), collapse=", ")))
    }
    if (!all(names(dt_new) == expected_cols)) {
      stop(paste("ERROR: Invalid columns in new file:", new_file, 
                 "\nExpected:", paste(expected_cols, collapse=", "), 
                 "\nFound:", paste(names(dt_new), collapse=", ")))
    }

    dt_old_filtered <- dt_old[population != "neural"]

    setorder(dt_old_filtered, population, anchor)
    setorder(dt_new, population, anchor)

    test_identity <- all.equal(dt_old_filtered, dt_new, check.attributes = FALSE)
    
    if (!isTRUE(test_identity)) {
      stop(paste0("ERROR: Data mismatch detected for time window ", tw, " and loop ", loop_id, "!\n",
                  "Details: ", paste(test_identity, collapse = "; ")))
    }
  }
  cat("Time window", tw, "verified successfully (all 417 loops match).\n")
}

cat("\n======================================================\n")
cat("SUCCESS: All files exist, have the correct columns,\n")
cat("and contain identical values (excluding 'neural') across pipelines!\n")
cat("======================================================\n")
