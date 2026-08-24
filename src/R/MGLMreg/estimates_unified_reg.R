library(data.table)

#####################################################################
# Command line arguments
#####################################################################

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop("Usage: Rscript extract_coefficients.R <time_set> <model_type> <init_type>\n",
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

output_file <- file.path(output_dir, paste0("coefficients_", time_set, "_", model_type, "_", init_type, ".tsv"))

loops_file <- file.path(datadir, "long_and_short_range_loops_D_mel.tsv")

if (!file.exists(loops_file)) {
  stop("Loops file not found: ", loops_file)
}

loops_table <- fread(loops_file)
loops_of_interest <- loops_table$loop_id
loops_of_interest <- loops_of_interest[order(as.integer(sub("^L", "", loops_of_interest)))]

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

  coef_mat <- fit_reg@coefficients
  if (is.null(coef_mat)) {
    warning("No coefficients matrix for loop ", loop_id)
    next
  }

  # Ensure row names correspond to design matrix predictor names
  if (!is.null(fit_reg@data$X) && ncol(fit_reg@data$X) == nrow(coef_mat)) {
    rownames(coef_mat) <- colnames(fit_reg@data$X)
  }

  dt <- as.data.table(coef_mat, keep.rownames = "parameter")
  dt[, loop_id := loop_id]

  setcolorder(dt, c("loop_id", "parameter", setdiff(names(dt), c("loop_id", "parameter"))))

  results_list[[loop_id]] <- dt
}

if (length(results_list) > 0) {
  results_dt <- rbindlist(results_list, fill = TRUE)
  fwrite(results_dt, output_file, sep = "\t")
  cat("Saved coefficients to:", output_file, "\n")
} else {
  cat("No coefficient results found for the selected loops.\n")
}
