library(dplyr)

config_path <- file.path("config", "config.yml")
config <- yaml::yaml.load_file(config_path)

datadir <- config$paths$datadir
resultsdir <- config$paths$resultsdir

input_file <- file.path(resultsdir, "lineage_analysis", "lineage_populations.tsv")
output_file <- file.path(resultsdir, "lineage_analysis", "time_window_clusters.txt")

# 1. Wczytanie wcześniej przygotowanego pliku TSV
df <- read.table(input_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)

# 2. Zdefiniowanie okien czasowych w ustalonej kolejności
time_windows <- c("00-02", "02-04", "04-06", "06-08", "08-10", 
                  "10-12", "12-14", "14-16", "16-18", "18-20")

# 3. Otwarcie połączenia do pliku tekstowego w trybie zapisu ("w" - write)
file_conn <- file(output_file, open = "w")

# 4. Pętla przechodząca przez każde okno czasowe
for (tw in time_windows) {
  
  # Zapisanie nagłówka okna czasowego (np. "00-02")
  writeLines(tw, file_conn)
  
  # Filtrowanie danych dla aktualnego okna czasowego
  df_subset <- df %>% filter(time == tw)
  
  # Jeśli w danym oknie czasowym są jakiekolwiek komórki/populacje
  if (nrow(df_subset) > 0) {
    # Stworzenie wektora z połączonymi napisami: "populacja, lineage"
    lines_to_write <- paste(df_subset$population, df_subset$lineage, sep = ", ")
    
    # Zapisanie utworzonych linii do pliku
    writeLines(lines_to_write, file_conn)
  }
}

# 5. Zamknięcie połączenia z plikiem
close(file_conn)

print(paste("Zestawienie klastrów i linii czasowych zostało zapisane w:", output_file))
