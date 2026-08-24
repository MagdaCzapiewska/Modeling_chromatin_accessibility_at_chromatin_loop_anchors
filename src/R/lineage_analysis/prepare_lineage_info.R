# Wczytanie niezbędnych bibliotek
library(readxl)
library(dplyr)
library(stringr)

config_path <- file.path("config", "config.yml")
config <- yaml::yaml.load_file(config_path)

datadir <- config$paths$datadir
resultsdir <- config$paths$resultsdir

input_file <- file.path(datadir, "calderon_data", "science.abn5800_table_s2.xlsx")

# ATAC cell type annotations and rationale. 
# ... [Twoje nagłówki i komentarze] ...

outdir <- file.path(resultsdir, "lineage_analysis")

# 1. Utworzenie katalogu wyjściowego, jeśli nie istnieje
if (!dir.exists(outdir)) {
  dir.create(outdir, recursive = TRUE)
}

# 2. Wczytanie pliku xlsx - pomijamy 1 wiersz, ponieważ nagłówki są w wierszu 2
df <- read_excel(input_file, skip = 1)

# 3. Przetwarzanie tabeli
df_processed <- df %>%
  # Wybór interesujących kolumn
  select(time, cluster, lineage, refined_annotation) %>%
  # Tworzenie nowej kolumny 'population'
  mutate(
    # Usunięcie kropek i zamiana spacji na "_"
    clean_annotation = str_replace_all(refined_annotation, "\\.", ""),
    clean_annotation = str_replace_all(clean_annotation, " ", "_"),
    
    # Sklejenie klastra z oczyszczoną adnotacją
    population = paste(cluster, clean_annotation, sep = "_")
  ) %>%
  # Pozostawienie tylko wymaganych kolumn
  select(time, cluster, lineage, refined_annotation, population)

# 4. Ścieżka do pliku wynikowego
output_file <- file.path(outdir, "lineage_populations.tsv")

# 5. Zapis do pliku .tsv
write.table(df_processed, file = output_file, sep = "\t", quote = FALSE, row.names = FALSE)

print(paste("Plik został pomyślnie zapisany w:", output_file))
