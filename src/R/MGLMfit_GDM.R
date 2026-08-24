library(data.table)
library(filelock)

config_path <- file.path("config", "config.yml")
config <- yaml::yaml.load_file(config_path)

srcdir <- config$paths$srcdir

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 6) {
  stop("Usage: Rscript MGLMfit_GDM.R <init> <mode> <pop_number> <input_file> <output_file_rds> <log_file>")
}

init_arg     <- args[1]  # default, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8
mode_arg     <- args[2]  # all / pop
pop_number   <- args[3]  # population number as string (e.g. "1", "2" or "NA" for synthetic data)
input_file   <- args[4]
rds_file     <- args[5]
log_file     <- args[6]

file_basename <- basename(input_file)

pkg_path <- file.path(srcdir, "MGLM")
r_files <- list.files(file.path(pkg_path, "R"), pattern = "\\.R$", full.names = TRUE)
for (f in r_files) source(f)

cat("Processing file:", file_basename, "| Mode:", mode_arg, "| Init:", init_arg, "\n")

if (!file.exists(input_file)) {
  log_line <- paste(
    format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    "ERROR_FILE_MISSING",
    paste("File:", file_basename),
    sep = " | "
  )
  lock <- lock(paste0(log_file, ".lock"))
  write(log_line, file = log_file, append = TRUE)
  unlock(lock)
  stop(log_line)
}

reads <- fread(
  input_file,
  select = c("x_A1", "x_A2", "total_reads", "population")
)

if (mode_arg == "pop") {
  if (!"population" %in% names(reads)) {
    stop("Column 'population' missing in input file, but mode='pop' was requested.")
  }

  reads[, pop_id := sub("^([0-9]+)_.*", "\\1", population)]

  reads <- reads[pop_id == pop_number]
  reads[, pop_id := NULL]
}

if ("population" %in% names(reads)) {
  reads[, population := NULL]
}

reads[, x_out := total_reads - x_A1 - x_A2]
reads[, total_reads := NULL]

y <- as.matrix(reads[, .(x_A1, x_A2, x_out)])

col_sums <- colSums(y)
ord <- order(col_sums, decreasing = TRUE)
y <- y[, ord, drop = FALSE]
col_sums <- col_sums[ord]

if (anyNA(reads$x_A1) || anyNA(reads$x_A2)) {
  log_line <- paste(
    format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    "SKIP_NA",
    paste("File:", file_basename),
    paste("Cells:", nrow(y)),
    "Reason: NA in counts",
    sep = " | "
  )

  lock <- lock(paste0(log_file, ".lock"))
  write(log_line, file = log_file, append = TRUE)
  unlock(lock)

  saveRDS(NULL, rds_file)
  quit(save = "no", status = 0)
}

col_sums <- colSums(y)
nonzero_cols <- sum(col_sums > 0)

rm(reads)
fit <- NULL

if (nonzero_cols < 3 || nrow(y) == 0) {
  zero_cols <- names(col_sums)[col_sums == 0]
  reason <- ifelse(nrow(y) == 0, "No cells", paste("ZeroCols:", paste(zero_cols, collapse = ",")))
  
  log_line <- paste(
    format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    "SKIP",
    paste("File:", file_basename),
    paste("Cells:", nrow(y)),
    paste("Reason:", reason),
    sep = " | "
  )
  
  lock <- lock(paste0(log_file, ".lock"))
  write(log_line, file = log_file, append = TRUE)
  unlock(lock)
  saveRDS(fit, rds_file)
  
} else {
  log_messages <- character()

  mglm_args <- list(data = y, dist = "GDM")
  if (init_arg != "default") {
    init_val <- as.numeric(init_arg)
    mglm_args$init <- matrix(rep(init_val, 4), nrow = 1, byrow = TRUE)
  }

  timing <- system.time(
    fit <- tryCatch({
      res <- withCallingHandlers({
        do.call(MGLMfit, mglm_args)
      }, warning = function(w) {
        log_messages <<- c(
          log_messages,
          paste(
            "WARNING",
            paste("File:", file_basename),
            "|",
            conditionMessage(w)
          )
        )
        invokeRestart("muffleWarning")
      })
      res
    }, error = function(e) {
      log_messages <<- c(
        log_messages,
        paste(
          "ERROR",
          paste("File:", file_basename),
          "|",
          e$message
        )
      )
      return(NULL)
    })
  )

  saveRDS(fit, rds_file)

  log_line <- paste(
    format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    if (is.null(fit)) "ERROR" else "SUCCESS",
    paste("File:", file_basename),
    paste("Cells:", nrow(y)),
    paste("ElapsedTime(s):", round(timing["elapsed"], 2)),
    sep = " | "
  )

  lock <- lock(paste0(log_file, ".lock"))
  write(log_line, file = log_file, append = TRUE)
  if (length(log_messages) > 0) {
    write(log_messages, file = log_file, append = TRUE)
  }
  unlock(lock)
}
