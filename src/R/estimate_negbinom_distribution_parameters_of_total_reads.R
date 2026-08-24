library(data.table)
library(MASS)
library(dplyr)
library(filelock)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  stop("Missing required arguments!\nUsage: Rscript estimate_negbinom_distribution_parameters_of_total_reads.R <tw>")
}

tw <- args[1]

config_path <- file.path("config", "config.yml")
config <- yaml::yaml.load_file(config_path)

datadir <- config$paths$datadir
resultsdir <- config$paths$resultsdir

total_reads_dir <- file.path(resultsdir, "counts", "total_reads")
output_dir <- file.path(resultsdir, "parameters", "negbinom_distribution_of_total_reads")

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
}

log_file <- file.path(output_dir, paste0("log_file_tw_", tw, ".txt"))

log_msg <- function(type = "INFO", tw, population, message_text, log_file) {
  msg <- paste0(
    format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    " | ", type, " | tw=", tw,
    " population=", population,
    " | ", message_text, "\n"
  )
  message(msg)

  lock_file_path <- paste0(log_file, ".lock")
  lock_handle <- filelock::lock(lock_file_path, timeout = Inf)
  
  cat(msg, file = log_file, append = TRUE)

  filelock::unlock(lock_handle)
}

safe_negbinom_fit <- function(reads_vector, population, tw, log_file) {
  
  output_df <- data.frame(
    population        = population,
    mu                = NA_real_,
    mu_se             = NA_real_,
    size              = NA_real_,
    size_se           = NA_real_,
    mu_est_over_se    = NA_real_,
    size_est_over_se  = NA_real_,
    min_est_over_se   = NA_real_,
    n_cells           = length(reads_vector)
  )
  
  fit_object <- NULL

  tryCatch({
    withCallingHandlers({
      
      fit_object <- MASS::fitdistr(reads_vector, densfun = "negative binomial")
      
      output_df$mu   <- fit_object$estimate[["mu"]]
      output_df$size <- fit_object$estimate[["size"]]

      output_df$mu_se   <- fit_object$sd[["mu"]]
      output_df$size_se <- fit_object$sd[["size"]]
      
      output_df$mu_est_over_se   <- output_df$mu / output_df$mu_se
      output_df$size_est_over_se <- output_df$size / output_df$size_se
      

      output_df$min_est_over_se  <- min(output_df$mu_est_over_se, output_df$size_est_over_se, na.rm = TRUE)
      
    }, warning = function(w) {
      log_msg("WARNING", tw, population, conditionMessage(w), log_file)
      invokeRestart("muffleWarning")
      
    }, error = function(e) {
      log_msg("ERROR", tw, population, conditionMessage(e), log_file)
    })
    
  }, error = function(e) {
    return(list(summary = output_df, fit = NULL))
  })
  
  return(list(summary = output_df, fit = fit_object))
}

input_reads_file <- file.path(total_reads_dir, paste0("total_reads_", tw, ".tsv.gz"))
if(!file.exists(input_reads_file)) {
  stop("Input file does not exist: ", input_reads_file)
}

dt <- fread(input_reads_file)

required_cols <- c("population", "total_reads")
missing_cols <- setdiff(required_cols, names(dt))
if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "), " in file: ", input_reads_file)
}

priors_list <- list()
models_list <- list()

############################################################################
# 1) ESTIMATION FOR GROUP: ALL
############################################################################
message("TW: ", tw, " - fitting ALL")

reads_all <- dt$total_reads

if (length(reads_all) == 0 || all(reads_all == 0)) {
  log_msg("INFO", tw, "all", "Skipping estimation: total_reads vector is empty or all elements are 0.", log_file)
  priors_list[[length(priors_list) + 1]] <- data.frame(
    population = "all", mu = NA_real_, mu_se = NA_real_, size = NA_real_, size_se = NA_real_,
    mu_est_over_se = NA_real_, size_est_over_se = NA_real_, min_est_over_se = NA_real_, n_cells = length(reads_all)
  )
  models_list[["all"]] <- NULL
} else {
  res_all <- safe_negbinom_fit(reads_all, "all", tw, log_file)
  priors_list[[length(priors_list) + 1]] <- res_all$summary
  models_list[["all"]] <- res_all$fit
}

############################################################################
# 2) ESTIMATION FOR EACH POPULATION SEPARATELY
############################################################################
pops <- unique(dt$population)

for (i in seq_along(pops)) {
  pop <- pops[i]
  message("TW: ", tw, " - fitting pop ", pop)
  
  reads_pop <- dt[population == pop, total_reads]
  
  if (length(reads_pop) == 0 || all(reads_pop == 0)) {
    log_msg("INFO", tw, pop, "Skipping estimation: total_reads vector is empty or all elements are 0.", log_file)
    priors_list[[length(priors_list) + 1]] <- data.frame(
      population = pop, mu = NA_real_, mu_se = NA_real_, size = NA_real_, size_se = NA_real_,
      mu_est_over_se = NA_real_, size_est_over_se = NA_real_, min_est_over_se = NA_real_, n_cells = length(reads_pop)
    )
    models_list[[pop]] <- NULL
  } else {
    res_pop <- safe_negbinom_fit(reads_pop, pop, tw, log_file)
    priors_list[[length(priors_list) + 1]] <- res_pop$summary
    models_list[[pop]] <- res_pop$fit
  }
}

############################################################################
# SAVE OUTPUTS (SPLIT INTO TSV.GZ AND RDS)
############################################################################

# Przygotowanie i posortowanie tabeli zbiorczej
final_priors <- bind_rows(priors_list) %>%
  mutate(cluster_id = ifelse(population == "all", -1, suppressWarnings(as.integer(sub("_.*", "", population))))) %>%
  mutate(cluster_id = ifelse(is.na(cluster_id), 999, cluster_id)) %>% 
  arrange(cluster_id) %>%
  select(-cluster_id)

# 1. Zapis tabeli statystyk do skompresowanego TSV
output_tsv_path <- file.path(output_dir, paste0("negbinom_parameters_", tw, "_all_and_by_population.tsv.gz"))
fwrite(final_priors, output_tsv_path, sep = "\t")

# 2. Zapis samych surowych obiektów modeli MASS do pliku RDS
output_rds_path <- file.path(output_dir, paste0("negbinom_models_", tw, "_all_and_by_population.rds"))
saveRDS(models_list, output_rds_path)

message("Job completed successfully for TW: ", tw)
