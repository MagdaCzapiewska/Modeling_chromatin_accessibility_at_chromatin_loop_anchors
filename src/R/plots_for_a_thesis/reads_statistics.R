#!/usr/bin/env Rscript

library(data.table)
library(yaml)

# 1. Wczytanie konfiguracji i ścieżek
config_path <- file.path("config", "config.yml")
if (!file.exists(config_path)) {
  stop("Błąd: Nie znaleziono pliku konfiguracji config/config.yml")
}

config <- yaml::yaml.load_file(config_path)
datadir <- config$paths$datadir
resultsdir <- config$paths$resultsdir

loops_file <- file.path(datadir, "long_and_short_range_loops_D_mel.tsv")
counts_dir <- file.path(resultsdir, "counts", "counts_in_anchors")

if (!file.exists(loops_file)) {
  stop("Błąd: Nie znaleziono pliku pętli: ", loops_file)
}

if (!dir.exists(counts_dir)) {
  stop("Błąd: Nie znaleziono katalogu z odczytami: ", counts_dir)
}

# 2. Pobranie listy pętli i okien czasowych
loops_table <- fread(loops_file)
if (!"loop_id" %in% names(loops_table)) {
  stop("Błąd: Kolumna 'loop_id' nie istnieje w pliku pętli.")
}
loop_ids <- loops_table$loop_id

time_windows <- c("00-02", "02-04", "04-06", "06-08", "08-10", 
                  "10-12", "12-14", "14-16", "16-18", "18-20")

# 3. Zmienne do akumulacji statystyk
total_cells <- 0
sum_total_reads <- 0
sum_x_A1 <- 0
sum_x_A2 <- 0
zeros_x_A1 <- 0
zeros_x_A2 <- 0
files_processed <- 0

message("=== Rozpoczynanie przetwarzania plików ===")

# 4. Pętla po oknach czasowych i pętlach
for (tw in time_windows) {
  for (loop_id in loop_ids) {
    
    # Sprawdzenie wersji skompresowanej i nieskompresowanej
    fpath <- file.path(counts_dir, paste0("reads_", tw, "_", loop_id, ".tsv.gz"))
    if (!file.exists(fpath)) {
      fpath <- file.path(counts_dir, paste0("reads_", tw, "_", loop_id, ".tsv"))
    }
    
    if (!file.exists(fpath)) next
    
    # Wczytujemy tylko potrzebne kolumny dla oszczędności pamięci i czasu
    dt <- fread(fpath, select = c("x_A1", "x_A2", "total_reads"))
    if (nrow(dt) == 0) next
    
    n <- nrow(dt)
    total_cells <- total_cells + n
    
    sum_total_reads <- sum_total_reads + sum(dt$total_reads, na.rm = TRUE)
    sum_x_A1        <- sum_x_A1 + sum(dt$x_A1, na.rm = TRUE)
    sum_x_A2        <- sum_x_A2 + sum(dt$x_A2, na.rm = TRUE)
    
    zeros_x_A1      <- zeros_x_A1 + sum(dt$x_A1 == 0, na.rm = TRUE)
    zeros_x_A2      <- zeros_x_A2 + sum(dt$x_A2 == 0, na.rm = TRUE)
    
    files_processed <- files_processed + 1
  }
}

# 5. Obliczenie statystyk końcowych
if (total_cells == 0) {
  stop("Nie znaleziono żadnych danych do przeanalizowania.")
}

mean_total_reads <- sum_total_reads / total_cells
mean_x_A1        <- sum_x_A1 / total_cells
mean_x_A2        <- sum_x_A2 / total_cells

pct_zeros_x_A1   <- (zeros_x_A1 / total_cells) * 100
pct_zeros_x_A2   <- (zeros_x_A2 / total_cells) * 100

# 6. Wypisanie wyników na konsolę
cat("\n==================================================\n")
cat("          PODSUMOWANIE STATYSTYK ODCZYTÓW         \n")
cat("==================================================\n")
cat(sprintf("Przetworzonych plików:           %s\n", format(files_processed, big.mark = ",")))
cat(sprintf("Łączna liczba obserwacji (N):    %s\n", format(total_cells, big.mark = ",")))
cat("--------------------------------------------------\n")
cat(sprintf("Średnia total_reads:            %.2f\n", mean_total_reads))
cat(sprintf("Średnia liczba odczytów x_A1:   %.2f\n", mean_x_A1))
cat(sprintf("Średnia liczba odczytów x_A2:   %.2f\n", mean_x_A2))
cat("--------------------------------------------------\n")
cat(sprintf("Procent zer w x_A1:              %.2f%%\n", pct_zeros_x_A1))
cat(sprintf("Procent zer w x_A2:              %.2f%%\n", pct_zeros_x_A2))
cat("==================================================\n\n")
