library(data.table)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage: Rscript naive_cor_real.R <tw> <mode> <output_file_tsv_gz> <real_counts_dir> [results_dir_for_pops]")
}

current_tw   <- args[1]
mode_arg     <- args[2] # "all" or "pops"
output_file  <- args[3]
counts_dir   <- args[4]

config <- yaml::yaml.load_file(file.path("config", "config.yml"))
loops_table <- fread(file.path(config$paths$datadir, "long_and_short_range_loops_D_mel.tsv"))
loop_ids <- unique(loops_table$loop_id)

cor_results <- list()

compute_naive <- function(reads_dt, meta_list) {
  base_dt <- as.data.table(meta_list)
  base_dt[, `:=`(
    spearman_rho = NA_real_, spearman_pvalue = NA_real_,
    pearson_r = NA_real_, pearson_pvalue = NA_real_
  )]
  
  if (nrow(reads_dt) < 3 || anyNA(reads_dt$x_A1) || anyNA(reads_dt$x_A2)) return(base_dt)
  
  ct_s <- suppressWarnings(tryCatch(cor.test(reads_dt$x_A1, reads_dt$x_A2, method = "spearman"), error = function(e) NULL))
  ct_p <- suppressWarnings(tryCatch(cor.test(reads_dt$x_A1, reads_dt$x_A2, method = "pearson"), error = function(e) NULL))
  
  if (!is.null(ct_s) && !is.null(ct_p)) {
    base_dt[, `:=`(
      spearman_rho = as.numeric(ct_s$estimate),
      spearman_pvalue = ct_s$p.value,
      pearson_r = as.numeric(ct_p$estimate),
      pearson_pvalue = ct_p$p.value
    )]
  }
  return(base_dt)
}

if (mode_arg == "all") {
  for (loop in loop_ids) {
    reads_file <- file.path(counts_dir, paste0("reads_", current_tw, "_", loop, ".tsv.gz"))
    if (!file.exists(reads_file)) next
    
    reads <- fread(reads_file, select = c("x_A1", "x_A2"))
    meta <- list(time_window = current_tw, loop_id = loop)
    cor_results[[length(cor_results) + 1]] <- compute_naive(reads, meta)
  }
} else {
  results_base_dir <- args[5]
  card_dt <- fread(file.path(results_base_dir, "cluster_cardinality.tsv.gz"))
  valid_pops <- as.character(card_dt[get(current_tw) > 0][[1]])
  
  for (loop in loop_ids) {
    reads_file <- file.path(counts_dir, paste0("reads_", current_tw, "_", loop, ".tsv.gz"))
    if (!file.exists(reads_file)) next
    
    reads_all <- fread(reads_file, select = c("x_A1", "x_A2", "population"))
    reads_all[, pop_id := sub("^([0-9]+)_.*", "\\1", population)]
    
    pop_mapping <- setNames(unique(reads_all$population), sub("^([0-9]+)_.*", "\\1", unique(reads_all$population)))
    
    for (pop in valid_pops) {
      sub_reads <- reads_all[pop_id == pop]
      full_pop_name <- if (pop %in% names(pop_mapping)) pop_mapping[pop] else paste0(pop, "_unknown")
      
      meta <- list(time_window = current_tw, population = full_pop_name, loop_id = loop)
      cor_results[[length(cor_results) + 1]] <- compute_naive(sub_reads[, .(x_A1, x_A2)], meta)
    }
  }
}

if (length(cor_results) > 0) {
  fwrite(rbindlist(cor_results, fill = TRUE), output_file, sep = "\t", compress = "gzip")
} else {
  empty <- data.table(time_window = current_tw, loop_id = NA_character_, spearman_rho = NA_real_)
  fwrite(empty, output_file, sep = "\t", compress = "gzip")
}
cat("Finished native correlation calculation for TW:", current_tw, "\n")
