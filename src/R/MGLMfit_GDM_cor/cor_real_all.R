library(data.table)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
  stop("Usage: Rscript cor_real_all.R <tw> <init> <output_file_tsv_gz> <fit_base_dir> <real_counts_dir>")
}

current_tw     <- args[1]
init_val       <- args[2]
output_file    <- args[3]
fit_base_dir   <- args[4]
counts_dir     <- args[5]

source(file.path("src", "R", "MGLMfit_GDM_cor", "correlation_functions.R"))

config <- yaml::yaml.load_file(file.path("config", "config.yml"))
loops_table <- fread(file.path(config$paths$datadir, "long_and_short_range_loops_D_mel.tsv"))
loop_ids <- unique(loops_table$loop_id)

fit_dir <- file.path(fit_base_dir, "MGLMfit_GDM", "real_data", "all", paste0("init_", init_val))
cor_results <- list()

for (loop in loop_ids) {
  reads_file <- file.path(counts_dir, paste0("reads_", current_tw, "_", loop, ".tsv.gz"))
  fit_file   <- file.path(fit_dir, paste0("fit_real_", current_tw, "_", loop, ".rds"))
  
  meta <- list(time_window = current_tw, loop_id = loop)
  
  if (!file.exists(reads_file)) next
  
  fit <- tryCatch(readRDS(fit_file), error = function(e) NULL)
  cor_results[[length(cor_results) + 1]] <- compute_gdm_correlations(fit, reads_file, mode = "all", metadata_list = meta)
}

if (length(cor_results) > 0) {
  fwrite(rbindlist(cor_results), output_file, sep = "\t", compress = "gzip")
} else {
  empty_template <- compute_gdm_correlations(NULL, NULL, metadata_list = list(time_window=current_tw, loop_id=NA_character_))
  fwrite(empty_template, output_file, sep = "\t", compress = "gzip")
}
