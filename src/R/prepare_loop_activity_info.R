library(data.table)
library(yaml)

# Zapobieganie zapisowi naukowemu (scipen) w całej sesji R
options(scipen = 999)

# 1. Wczytanie konfiguracji i ścieżek
config_path <- file.path("config", "config.yml")
if (!file.exists(config_path)) {
  stop("Nie znaleziono pliku konfiguracji: config/config.yml")
}

config <- yaml::yaml.load_file(config_path)

datadir    <- config$paths$datadir
srcdir     <- config$paths$srcdir
resultsdir <- config$paths$resultsdir

dir.create(resultsdir, showWarnings = FALSE, recursive = TRUE)

loops_file  <- file.path(datadir, "long_and_short_range_loops_D_mel.tsv")
output_file <- file.path(resultsdir, "loop_activity_info.tsv")

if (!file.exists(loops_file)) {
  stop("Nie znaleziono pliku z pętlami: ", loops_file)
}

# 2. Wczytanie danych
loops_data <- fread(loops_file)

# 3. Definicja kolumn aktywności
activity_cols <- c(
  "Dmel_6-8h_Neuroblasts",
  "Dmel_6-8h_Neurons",
  "Dmel_6-8h_Glia",
  "Dmel_10-12h_Neuroblasts",
  "Dmel_10-12h_Neurons",
  "Dmel_10-12h_Glia",
  "Dmel_14-16h_Neuroblasts",
  "Dmel_14-16h_Neurons",
  "Dmel_14-16h_Glia",
  "Dmel_16-18h_whole_embryo",
  "Dmel_larval_brain",
  "Dmel_adult_brain"
)

# Sprawdzenie obecności kolumn w pliku wejściowym
missing_cols <- setdiff(c("loop_id", activity_cols), names(loops_data))
if (length(missing_cols) > 0) {
  stop("Brakujące kolumny w pliku wejściowym: ", paste(missing_cols, collapse = ", "))
}

# 4. Wyciągnięcie wybranych kolumn
res_dt <- loops_data[, c("loop_id", activity_cols), with = FALSE]

# 5. Bezpośrednia konwersja kolumn na typ tekstowy (character) i zamiana NA na "n"
char_dt <- res_dt[, lapply(.SD, function(x) fifelse(is.na(x), "n", as.character(x))), .SDcols = activity_cols]

# 6. Sklejenie w profile
res_dt[, small_profile := do.call(paste0, char_dt[, 1:9, with = FALSE])]
res_dt[, large_profile := do.call(paste0, char_dt[, 1:12, with = FALSE])]

# 7. Zapis z wymuszeniem cudzysłowów (quote = TRUE)
fwrite(res_dt, output_file, sep = "\t", quote = TRUE)
cat("[INFO] Zapisano plik z profilami aktywności do:", output_file, "\n")
