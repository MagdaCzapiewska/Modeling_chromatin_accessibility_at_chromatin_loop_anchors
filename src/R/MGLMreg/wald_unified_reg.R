library(data.table)

#####################################################################
# Command line arguments
#####################################################################

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop("Usage: Rscript wald_unified_reg.R <time_set> <model_type> <init_type>\n",
       "  time_set:   '0h+' or '10h+'\n",
       "  model_type: 'time' or 'time_tissue'\n",
       "  init_type:  'smart' or 'default'")
}

time_set   <- args[1]  # "0h+" | "10h+"
model_type <- args[2]  # "time" | "time_tissue"
init_type  <- args[3]  # "smart" | "default"

cat("Time set:", time_set, "\n")
cat("Model type:", model_type, "\n")
cat("Init type:", init_type, "\n")

# Path configuration
config_path <- file.path("config", "config.yml")
config <- yaml::yaml.load_file(config_path)

datadir    <- config$paths$datadir
resultsdir <- config$paths$resultsdir

input_dir   <- file.path(resultsdir, "MGLMreg", paste0("GDM_", time_set, "_", model_type, "_", init_type))
output_dir  <- file.path(resultsdir, "MGLMreg")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

output_file <- file.path(output_dir, paste0("wald_test_results_", time_set, "_", model_type, "_", init_type, ".tsv"))

loops_file <- file.path(datadir, "long_and_short_range_loops_D_mel.tsv")

if (!file.exists(loops_file)) {
  stop("Loops file not found: ", loops_file)
}

loops_table <- fread(loops_file)
loops_of_interest <- loops_table$loop_id
loops_of_interest <- loops_of_interest[order(as.integer(sub("^L", "", loops_of_interest)))]

# Helper function to extract values safely
get_coef_safe <- function(mat, row, col) {
  if (!is.null(mat) && row %in% rownames(mat) && col %in% colnames(mat)) {
    return(mat[row, col])
  } else {
    return(NA_real_)
  }
}

results_list <- list()

for (loop_id in loops_of_interest) {
  cat("Processing loop:", loop_id, "\n")
  file <- file.path(input_dir, paste0("fit_reg_", loop_id, ".rds"))

  if (!file.exists(file)) {
    warning("Missing file for loop ", loop_id)
    next
  }

  fit_reg <- readRDS(file)

  if (is.null(fit_reg)) {
    warning("fit_reg is NULL for loop ", loop_id)
    next
  }

  test_mat <- fit_reg@test
  if (is.null(test_mat)) {
    warning("No Wald test results for loop ", loop_id)
    next
  }

  # Obliczanie statusu dopasowania (fit_status)
  grad_mat <- fit_reg@gradient
  is_converged <- !is.null(grad_mat) && (mean(grad_mat^2, na.rm = TRUE) <= 1e-04)
  p_val_time <- get_coef_safe(test_mat, "time", "Pr(>wald)")

  if (is.na(p_val_time)) {
    fit_status <- "FAILED_PVAL_NA"
  } else if (!is_converged) {
    fit_status <- "NOT_CONVERGED"
  } else {
    fit_status <- "SUCCESS"
  }

  dt <- as.data.table(test_mat, keep.rownames = "parameter")
  dt[, loop_id := loop_id]
  dt[, fit_status := fit_status]

  setcolorder(dt, c("loop_id", "parameter", "fit_status", setdiff(names(dt), c("loop_id", "parameter", "fit_status"))))

  results_list[[loop_id]] <- dt
}

if (length(results_list) > 0) {
  results_dt <- rbindlist(results_list, fill = TRUE)
  alpha_threshold <- 0.05
  
  # Poprawka BH aplikowana WYŁĄCZNIE dla udanych dopasowań (SUCCESS)
  #if ("Pr(>wald)" %in% names(results_dt)) {
  #  results_dt[fit_status == "SUCCESS" & !is.na(`Pr(>wald)`), 
  #             padj_BH := p.adjust(`Pr(>wald)`, method = "BH"), 
  #             by = parameter]
  #  results_dt[, is_significant := !is.na(padj_BH) & padj_BH < alpha_threshold]
  #}
  # Poprawka BH aplikowana globalnie dla WSZYSTKICH udanych dopasowań (SUCCESS)
  if ("Pr(>wald)" %in% names(results_dt)) {
    results_dt[fit_status == "SUCCESS" & !is.na(`Pr(>wald)`), 
               padj_BH := p.adjust(`Pr(>wald)`, method = "BH")]
    
    results_dt[, is_significant := !is.na(padj_BH) & padj_BH < alpha_threshold]
  }

  fwrite(results_dt, output_file, sep = "\t")
  cat("Saved Wald test results with fit status and filtered BH correction to:", output_file, "\n")
} else {
  cat("No Wald test results found for the selected loops.\n")
}
