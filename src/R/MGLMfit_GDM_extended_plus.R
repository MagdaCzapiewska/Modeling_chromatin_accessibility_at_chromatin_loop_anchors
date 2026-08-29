library(data.table)
library(filelock)

config_path <- file.path("config", "config.yml")
config <- yaml::yaml.load_file(config_path)
srcdir <- config$paths$srcdir

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 6) {
  stop("Usage: Rscript MGLMfit_GDM_extended_plus.R <init> <rho> <seed> <synthetic_dir> <output_base_dir> <log_file>")
}

init_arg        <- args[1]
rho_val         <- as.numeric(args[2])
seed_val        <- as.integer(args[3])
synthetic_dir   <- args[4]
output_base_dir <- args[5]
log_file        <- args[6]

rho_str <- format(rho_val, nsmall = 1)

# Nowa rozszerzona siatka parametrów
N_CELLS        <- c(1000, 2000, 5000, 10000, 20000)
MUS            <- c(1000, 2000, 3000, 4000, 5000)
SIZE_NEGBINOMS <- c("0.1", "0.2", "0.5", "1.0", "2.0", "inf", "fixed")

ALPHA <- 0.15
BETA  <- 600

pkg_path <- file.path(srcdir, "MGLM")
r_files <- list.files(file.path(pkg_path, "R"), pattern = "\\.R$", full.names = TRUE)
for (f in r_files) source(f)

target_out_dir <- file.path(output_base_dir, paste0("rho_", rho_str), paste0("init_", init_arg))
if (!dir.exists(target_out_dir)) {
  dir.create(target_out_dir, recursive = TRUE, showWarnings = FALSE)
}

# Ścieżka do jednego zbiorczego pliku TSV dla danego ziarna
input_file <- file.path(synthetic_dir, paste0("rho_", rho_str), paste0("synthetic_counts_seed", seed_val, ".tsv.gz"))

if (!file.exists(input_file)) {
  log_line <- paste(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "ERROR_FILE_MISSING", paste("File:", input_file), sep = " | ")
  lock <- lock(paste0(log_file, ".lock"))
  write(log_line, file = log_file, append = TRUE)
  unlock(lock)
  stop(paste("Input file missing:", input_file))
}

# Wczytanie całego pliku dla danego seeda
all_reads <- fread(input_file)
setkey(all_reads, n_cells, mu, size_nb)

for (n in N_CELLS) {
  for (mu_val in MUS) {
    for (size_nb_str in SIZE_NEGBINOMS) {
      
      file_tag <- paste0("n", n, "_mu", mu_val, "_sizeNB", size_nb_str, "_alpha", ALPHA, "_beta", BETA, "_seed", seed_val)
      rds_file <- file.path(target_out_dir, paste0("fit_synthetic_", file_tag, ".rds"))
      
      # Wyciągnięcie konkretnego scenariusza z wygenerowanej ramki
      #reads <- all_reads[n_cells == n & mu == mu & size_nb == size_nb, .(x_A1, x_A2, total_reads)]
      #reads <- all_reads[n_cells == n & mu == ..mu_val & size_nb == ..size_nb_str, .(x_A1, x_A2, total_reads)]
      reads <- all_reads[.(n, mu_val, size_nb_str), .(x_A1, x_A2, total_reads), nomatch = NULL]
      
      if (nrow(reads) == 0) {
        log_line <- paste(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "SKIP_NO_DATA", paste("Scenario:", file_tag), sep = " | ")
        lock <- lock(paste0(log_file, ".lock"))
        write(log_line, file = log_file, append = TRUE)
        unlock(lock)
        saveRDS(NULL, rds_file)
        next
      }
      
      reads[, x_out := total_reads - x_A1 - x_A2]
      reads[, total_reads := NULL]
      
      y <- as.matrix(reads[, .(x_A1, x_A2, x_out)])
      
      col_sums <- colSums(y)
      ord <- order(col_sums, decreasing = TRUE)
      y <- y[, ord, drop = FALSE]
      col_sums <- col_sums[ord]
      
      if (anyNA(y) || any(y < 0)) {
        log_line <- paste(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "SKIP_INVALID", paste("Scenario:", file_tag), paste("Cells:", nrow(y)), "Reason: NA or negative counts in y", sep = " | ")
        lock <- lock(paste0(log_file, ".lock"))
        write(log_line, file = log_file, append = TRUE)
        unlock(lock)
        saveRDS(NULL, rds_file)
        next
      }
      
      nonzero_cols <- sum(col_sums > 0)
      rm(reads)
      fit <- NULL
      
      if (nonzero_cols < 3 || nrow(y) == 0) {
        zero_cols <- names(col_sums)[col_sums == 0]
        reason <- ifelse(nrow(y) == 0, "No cells", paste("ZeroCols:", paste(zero_cols, collapse = ",")))
        
        log_line <- paste(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "SKIP", paste("Scenario:", file_tag), paste("Cells:", nrow(y)), paste("Reason:", reason), sep = " | ")
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
              log_messages <<- c(log_messages, paste("WARNING", paste("Scenario:", file_tag), "|", conditionMessage(w)))
              invokeRestart("muffleWarning")
            })
            res
          }, error = function(e) {
            log_messages <<- c(log_messages, paste("ERROR", paste("Scenario:", file_tag), "|", e$message))
            return(NULL)
          })
        )
        
        saveRDS(fit, rds_file)
        
        log_line <- paste(
          format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
          if (is.null(fit)) "ERROR" else "SUCCESS",
          paste("Scenario:", file_tag),
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
    }
  }
}

done_file <- file.path(target_out_dir, paste0("fit_seed", seed_val, ".done"))
file.create(done_file)
cat("Successfully completed MGLM GDM fitting for rho:", rho_str, "init:", init_arg, "seed:", seed_val, "\n")
