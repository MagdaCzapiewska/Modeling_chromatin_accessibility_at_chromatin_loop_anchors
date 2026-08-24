library(data.table)
library(MGLM)
library(filelock)

#####################################################################

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 4) {
  stop("Usage: Rscript MN_time_0h+.R <loop_id> <out_file> <all> <n_take>")
}

loop_id <- args[1]
out_file <- args[2]
all_flag <- as.logical(args[3])
n_take <- as.numeric(args[4])

cat("Loop ID: ", loop_id, "\n")
cat("Output file: ", out_file, "\n")
cat("Use all reads:", all_flag, "\n")
cat("n_take:", n_take, "\n")

config_path <- file.path("config", "config.yml")
config <- yaml::yaml.load_file(config_path)

datadir <- config$paths$datadir
resultsdir <- config$paths$resultsdir

input_dir <- file.path(resultsdir, "counts", "counts_in_anchors")
output_dir <- dirname(out_file)
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

log_file <- file.path(output_dir, "MN_time_0h+_reg_all.log")

time_windows <- c("00-02", "02-04", "04-06", "06-08", "08-10", "10-12", "12-14", "14-16", "16-18", "18-20")

df_list <- list()

for (tw in time_windows) {

  f <- file.path(input_dir, paste0("reads_", tw, "_", loop_id, ".tsv.gz"))
  if (!file.exists(f)) {
    cat("Missing file for window:", tw, "- skip\n")
    next
  }

  reads <- fread(
    f,
    select = c("x_A1", "x_A2", "total_reads", "NNv1_age")
  )

  if (!all_flag && nrow(reads) > n_take) {
    set.seed(123)
    reads <- reads[sample(.N, n_take)]
  }

  reads[, x_out := total_reads - x_A1 - x_A2]

  df_list[[tw]] <- reads[, .(
    x_A1,
    x_A2,
    x_out,
    time = NNv1_age
  )]

  rm(reads)
  #gc()
}

df <- rbindlist(df_list)

if (nrow(df) == 0) {
  log_line <- paste(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "| SKIP | No data for loop", loop_id)
  
  lock <- lock(paste0(log_file, ".lock"))
  write(log_line, file = log_file, append = TRUE)
  unlock(lock)
  
  saveRDS(NULL, out_file)
  quit(save = "no", status = 0)
}

log_messages <- character()
fit_reg <- NULL

time_fit <- system.time({
  fit_reg <- tryCatch({
    
    res <- withCallingHandlers({
      
      MGLMreg(cbind(x_out, x_A2, x_A1) ~ time, data = df, dist = "MN")
      
    }, warning = function(w) {
      log_messages <<- c(
        log_messages,
        paste(
          format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
          "| WARNING | Loop:", loop_id,
          "|", conditionMessage(w)
        )
      )
      invokeRestart("muffleWarning")
    })
    
    res
    
  }, error = function(e) {
    log_messages <<- c(
      log_messages,
      paste(
        format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
        "| ERROR | Loop:", loop_id,
        "|", e$message
      )
    )
    return(NULL)
  })
})


status_tag <- if (is.null(fit_reg)) "ERROR_FIT_FAILED" else "SUCCESS"
log_line_main <- paste(
  format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
  status_tag,
  paste("Loop:", loop_id),
  paste("Rows:", nrow(df)),
  paste("ElapsedTime(s):", round(time_fit["elapsed"], 2)),
  sep = " | "
)

lock <- lock(paste0(log_file, ".lock"))

write(log_line_main, file = log_file, append = TRUE)

if (length(log_messages) > 0) {
  write(log_messages, file = log_file, append = TRUE)
}

unlock(lock)

if (!is.null(fit_reg)) {
  saveRDS(fit_reg, out_file)
  cat(paste(Sys.time(), "Loop", loop_id, "fitting succeeded. Model saved to:", out_file, "\n"))
  print(fit_reg)
} else {
  saveRDS(NULL, out_file)
  cat(paste(Sys.time(), "Loop", loop_id, "fitting FAILED.\n"))
}
