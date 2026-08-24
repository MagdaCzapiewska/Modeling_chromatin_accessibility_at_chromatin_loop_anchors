library(data.table)
library(ebbr)
library(dplyr)
library(filelock)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("Missing required arguments!\nUsage: Rscript estimate_beta_distribution_parameters_of_counts_in_anchors.R <tw> <loop_id>")
}

tw      <- args[1]
loop_id <- args[2]

config_path <- file.path("config", "config.yml")
config <- yaml::yaml.load_file(config_path)

datadir <- config$paths$datadir
resultsdir <- config$paths$resultsdir

counts_dir <- file.path(resultsdir, "counts", "counts_in_anchors")
output_dir <- file.path(resultsdir, "parameters", "beta_distribution_of_counts_in_anchors")

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
}

log_file <- file.path(output_dir, paste0("log_file_tw_", tw, ".txt"))

log_msg <- function(type = "INFO", loop_id, anchor, tw, population, message_text, log_file) {
  msg <- paste0(
    format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    " | ", type, " | loop_id=", loop_id,
    " anchor=", anchor,
    " tw=", tw,
    " population=", population,
    " | ", message_text, "\n"
  )
  message(msg)

  lock_file_path <- paste0(log_file, ".lock")

  lock_handle <- filelock::lock(lock_file_path, timeout = Inf)
  
  cat(msg, file = log_file, append = TRUE)

  filelock::unlock(lock_handle)
}

safe_ebb_fit_prior <- function(ebb_data, population, anchor, tw, loop_id, log_file) {
  output_df <- data.frame(
    population   = population,
    anchor       = anchor,
    alpha_prior  = NA_real_,
    beta_prior   = NA_real_,
    n_cells      = nrow(ebb_data)
  )

  tryCatch({
    withCallingHandlers({
      fit <- ebb_fit_prior(tbl = ebb_data, x = successes, n = trials, method = "mle")
      
      output_df$alpha_prior <- fit$parameters$alpha
      output_df$beta_prior  <- fit$parameters$beta
      
    }, warning = function(w) {
      log_msg("WARNING", loop_id, anchor, tw, population, conditionMessage(w), log_file)
      invokeRestart("muffleWarning") # Suppress warning propagation after logging
      
    }, error = function(e) {
      log_msg("ERROR", loop_id, anchor, tw, population, conditionMessage(e), log_file)
      # Error will automatically fall through to tryCatch handler
    })
    
  }, error = function(e) {
    # If a critical error happened, output_df retains NA values
    return(output_df)
  })
  
  return(output_df)
}

input_reads_file <- file.path(counts_dir, paste0("reads_", tw, "_", loop_id, ".tsv.gz"))
if(!file.exists(input_reads_file)) {
  stop("Input file does not exist: ", input_reads_file)
}

dt <- fread(input_reads_file)

loop_cols <- c("x_A1", "x_A2")
present_cols <- intersect(names(dt), loop_cols)
if (length(present_cols) == 0) {
  stop("Columns x_A1 or x_A2 are missing in file: ", input_reads_file)
}

priors_list <- list()

############################################################################
# 1) PRIOR ESTIMATION FOR GROUP: ALL
############################################################################
message("TW: ", tw, ", LOOP: ", loop_id,  " - fitting ALL")

for (anchor in c("A1", "A2")) {
  colname <- paste0("x_", anchor)
  if (!(colname %in% names(dt))) next
  
  ebb_data <- data.frame(successes = dt[[colname]], trials = dt$total_reads)
  
  if (all(ebb_data$successes == 0) || nrow(ebb_data) == 0) {
    log_msg("INFO", loop_id, anchor, tw, "all", "Skipping prior estimation: all successes are 0 or data is empty.", log_file)
    priors_list[[length(priors_list) + 1]] <- data.frame(
      population="all", anchor=anchor,
      alpha_prior=NA_real_, beta_prior=NA_real_,
      n_cells=nrow(ebb_data)
    )
  } else {
    priors_list[[length(priors_list) + 1]] <- safe_ebb_fit_prior(ebb_data, "all", anchor, tw, loop_id, log_file)
  }
}

############################################################################
# 2) PRIOR ESTIMATION FOR EACH POPULATION SEPARATELY
############################################################################
pops <- unique(dt$population)

for (i in seq_along(pops)) {
  pop <- pops[i]
  message("TW: ", tw, ", LOOP: ", loop_id,  " - fitting pop ", pop)
  dt_pop <- dt[population == pop]
  
  for (anchor in c("A1", "A2")) {
    colname <- paste0("x_", anchor)
    if (!(colname %in% names(dt_pop))) next
    
    ebb_data <- data.frame(successes = dt_pop[[colname]], trials = dt_pop$total_reads)
    
    if (all(ebb_data$successes == 0) || nrow(ebb_data) == 0) {
      log_msg("INFO", loop_id, anchor, tw, pop, "Skipping prior estimation: all successes are 0 or data is empty.", log_file)
      priors_list[[length(priors_list) + 1]] <- data.frame(
        population=pop, anchor=anchor,
        alpha_prior=NA_real_, beta_prior=NA_real_,
        n_cells=nrow(ebb_data)
      )
    } else {
      priors_list[[length(priors_list) + 1]] <- safe_ebb_fit_prior(ebb_data, pop, anchor, tw, loop_id, log_file)
    }
  }
}

############################################################################
# SAVE OUTPUTS
############################################################################
final_priors <- bind_rows(priors_list) %>%
  mutate(cluster_id = ifelse(population == "all", -1, as.integer(sub("_.*", "", population)))) %>%
  arrange(cluster_id) %>%
  select(-cluster_id)

output_file_path <- file.path(output_dir, paste0("beta_parameters_", tw, "_", loop_id, "_all_and_by_population.tsv.gz"))
fwrite(final_priors, output_file_path, sep="\t")
