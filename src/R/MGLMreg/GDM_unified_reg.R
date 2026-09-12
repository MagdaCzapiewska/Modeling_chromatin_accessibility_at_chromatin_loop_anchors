library(data.table)
library(MGLM)
library(filelock)

#####################################################################
# Command line arguments
#####################################################################

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 7) {
  stop("Usage: Rscript GDM_unified_reg.R <loop_id> <out_file> <all> <n_take> <time_set> <model_type> <init_type>\n",
       "  time_set:   '0h+' or '10h+'\n",
       "  model_type: 'time' or 'time_tissue'\n",
       "  init_type:  'smart' or 'default'")
}

loop_id    <- args[1]
out_file   <- args[2]
all_flag   <- as.logical(args[3])
n_take     <- as.numeric(args[4])
time_set   <- args[5]  # "0h+" | "10h+"
model_type <- args[6]  # "time" | "time_tissue"
init_type  <- args[7]  # "smart" | "default"

cat("Loop ID:", loop_id, "\n")
cat("Output file:", out_file, "\n")
cat("Use all reads:", all_flag, "\n")
cat("n_take:", n_take, "\n")
cat("Time set:", time_set, "\n")
cat("Model type:", model_type, "\n")
cat("Init type:", init_type, "\n")

# Path configuration
config_path <- file.path("config", "config.yml")
config <- yaml::yaml.load_file(config_path)

resultsdir <- config$paths$resultsdir
srcdir     <- config$paths$srcdir

input_dir  <- file.path(resultsdir, "counts", "counts_in_anchors")
output_dir <- dirname(out_file)
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

log_file <- file.path(output_dir, paste0("GDM_", time_set, "_", model_type, "_", init_type, ".log"))

# Load local MGLM package source if available
pkg_path <- file.path(srcdir, "MGLM")
if (dir.exists(pkg_path)) {
  r_files <- list.files(file.path(pkg_path, "R"), pattern = "\\.R$", full.names = TRUE)
  for (f in r_files) source(f)
}

# Select time windows
if (time_set == "0h+") {
  time_windows <- c("00-02", "02-04", "04-06", "06-08", "08-10", "10-12", "12-14", "14-16", "16-18", "18-20")
} else if (time_set == "10h+") {
  time_windows <- c("10-12", "12-14", "14-16", "16-18", "18-20")
} else {
  stop("Invalid time_set argument. Choose '0h+' or '10h+'.")
}

# Select columns to read
select_cols <- c("x_A1", "x_A2", "total_reads", "NNv1_age")
if (model_type == "time_tissue") {
  select_cols <- c(select_cols, "population")
}

# Read data
df_list <- list()

for (tw in time_windows) {
  f <- file.path(input_dir, paste0("reads_", tw, "_", loop_id, ".tsv.gz"))
  if (!file.exists(f)) {
    cat("Missing file for window:", tw, "- skip\n")
  }

  reads <- fread(f, select = select_cols)

  if (!all_flag && nrow(reads) > n_take) {
    set.seed(123)
    reads <- reads[sample(.N, n_take)]
  }

  reads[, x_out := total_reads - x_A1 - x_A2]

  if (model_type == "time_tissue") {
    reads[, tissue := sub("^[^_]+_", "", population)]
    reads[, tissue := factor(tissue, levels = c("Unknown", setdiff(unique(tissue), "Unknown")))]
    reads[, tissue := droplevels(tissue)]

    df_list[[tw]] <- reads[, .(
      x_A1, x_A2, x_out,
      time = NNv1_age,
      tissue
    )]
  } else {
    df_list[[tw]] <- reads[, .(
      x_A1, x_A2, x_out,
      time = NNv1_age
    )]
  }

  rm(reads)
}

if (length(df_list) == 0) {
  log_line <- paste(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "| SKIP | No data for loop", loop_id)
  lock <- lock(paste0(log_file, ".lock"))
  write(log_line, file = log_file, append = TRUE)
  unlock(lock)
  saveRDS(NULL, out_file)
  quit(save = "no", status = 0)
}

df <- rbindlist(df_list)

if (nrow(df) == 0) {
  log_line <- paste(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "| SKIP | Empty table for loop", loop_id)
  lock <- lock(paste0(log_file, ".lock"))
  write(log_line, file = log_file, append = TRUE)
  unlock(lock)
  saveRDS(NULL, out_file)
  quit(save = "no", status = 0)
}

# Handle tissue factor levels
if (model_type == "time_tissue") {
  df[, tissue := droplevels(tissue)]
  if ("Unknown" %in% levels(df$tissue)) {
    df[, tissue := relevel(tissue, ref = "Unknown")]
  } else {
    warning("Unknown not present in data — baseline changed automatically.")
  }
}

# Define model formula and initialization matrix
if (model_type == "time") {
  fit_formula <- cbind(x_out, x_A2, x_A1) ~ time
  rhs_formula <- ~ time
} else if (model_type == "time_tissue") {
  fit_formula <- cbind(x_out, x_A2, x_A1) ~ time + tissue
  rhs_formula <- ~ time + tissue
} else {
  stop("Invalid model_type argument. Choose 'time' or 'time_tissue'.")
}

if (init_type == "smart") {
Y_temp <- as.matrix(df[, .(x_out, x_A2, x_A1)])
  X_temp <- model.matrix(rhs_formula, data = df)
  n_pred <- ncol(X_temp)
  n_cat  <- ncol(Y_temp)

  init_values <- matrix(0, nrow = n_pred, ncol = 2 * (n_cat - 1))

  init_values[1, ] <- c(9, 0.5, -0.5, -0.5)
} else if (init_type == "default") {
  init_values <- NULL
} else {
  stop("Invalid init_type argument. Choose 'smart' or 'default'.")
}

# Model fitting
log_messages <- character()
fit_reg <- NULL

time_fit <- system.time({
  fit_reg <- tryCatch({
    res <- withCallingHandlers({
      if (is.null(init_values)) {
        MGLMreg(fit_formula, data = df, dist = "GDM")
      } else {
        MGLMreg(fit_formula, data = df, dist = "GDM", init = init_values)
      }
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

# Save results and log output
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
