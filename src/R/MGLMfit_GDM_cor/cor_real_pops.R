library(data.table)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
  stop("Usage: Rscript cor_real_pops.R <tw> <init> <output_file_tsv_gz> <fit_base_dir> <real_counts_dir>")
}

current_tw     <- args[1]
init_val       <- args[2]
output_file    <- args[3]
fit_base_dir   <- args[4]
counts_dir     <- args[5]

source(file.path("src", "R", "MGLMfit_GDM_cor", "correlation_functions.R"))

config <- yaml::yaml.load_file(file.path("config", "config.yml"))

cardinality_file <- file.path(fit_base_dir, "cluster_cardinality.tsv.gz")
card_dt <- fread(cardinality_file)

if (!current_tw %in% names(card_dt)) {
  stop(paste("Error: Time window", current_tw, "not found in cluster_cardinality columns!"))
}

valid_pops <- as.character(card_dt[get(current_tw) > 0][[1]])

loops_table <- fread(file.path(config$paths$datadir, "long_and_short_range_loops_D_mel.tsv"))
loop_ids <- unique(loops_table$loop_id)

fit_dir <- file.path(fit_base_dir, "MGLMfit_GDM", "real_data", "pops", paste0("init_", init_val))
cor_results <- list()

for (loop in loop_ids) {
  reads_file <- file.path(counts_dir, paste0("reads_", current_tw, "_", loop, ".tsv.gz"))
  if (!file.exists(reads_file)) next

  pop_mapping <- tryCatch({
    pop_col <- fread(reads_file, select = "population")
    unique_pops <- unique(pop_col$population)
    keys <- sub("^([0-9]+)_.*", "\\1", unique_pops)
    setNames(unique_pops, keys)
  }, error = function(e) {
    NULL
  })

  for (pop in valid_pops) {
    fit_file <- file.path(fit_dir, paste0("fit_real_", current_tw, "_", loop, "_pop_", pop, ".rds"))

    full_pop_name <- if (!is.null(pop_mapping) && pop %in% names(pop_mapping)) {
      pop_mapping[pop]
    } else {
      paste0(pop, "_unknown")
    }
    
    meta <- list(time_window = current_tw, population = full_pop_name, loop_id = loop)
    
    fit <- tryCatch(readRDS(fit_file), error = function(e) NULL)
    cor_results[[length(cor_results) + 1]] <- compute_gdm_correlations(fit, reads_file, mode = "pop", pop_number = pop, metadata_list = meta)
  }
}

if (length(cor_results) > 0) {
  fwrite(rbindlist(cor_results), output_file, sep = "\t", compress = "gzip")
} else {
  empty_template <- compute_gdm_correlations(NULL, NULL, mode = "pop", metadata_list = list(time_window=current_tw, population=NA_character_, loop_id=NA_character_))
  fwrite(empty_template, output_file, sep = "\t", compress = "gzip")
}
