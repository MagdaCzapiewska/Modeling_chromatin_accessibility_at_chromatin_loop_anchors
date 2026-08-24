# Wczytanie niezbędnych bibliotek
library(dplyr)
library(stringr)

# Wczytanie ścieżek z konfiguracji YAML
config_path <- file.path("config", "config.yml")
config <- yaml::yaml.load_file(config_path)

resultsdir <- config$paths$resultsdir
input_file <- file.path(resultsdir, "lineage_analysis", "lineage_populations.tsv")

# Zdefiniowanie ścieżek dla wszystkich plików wynikowych
output_file_full    <- file.path(resultsdir, "lineage_analysis", "cluster_tree_indented.txt")
output_file_lineage <- file.path(resultsdir, "lineage_analysis", "lineage_only_tree_indented.txt")
output_file_parents <- file.path(resultsdir, "lineage_analysis", "lineage_parents_mapping.tsv")
output_file_paths   <- file.path(resultsdir, "lineage_analysis", "root_to_leaf_paths.txt")

# 1. Wczytanie danych wejściowych i przygotowanie struktury
df <- read.table(input_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
df$row_id <- 1:nrow(df)
df$parent_id <- NA_integer_  # Kolumna przechowująca relacje nadrzędne

# Chronologiczne okna czasowe
time_windows <- c("00-02", "02-04", "04-06", "06-08", "08-10", 
                  "10-12", "12-14", "14-16", "16-18", "18-20")

# 2. ALGORYTM: Szukanie rodzica z inteligentną logiką wyboru rodzeństwa (alfabetyczny fallback)
for (i in 2:length(time_windows)) {
  current_tw <- time_windows[i]
  previous_tw <- time_windows[i-1]
  
  current_rows <- df[df$time == current_tw, ]
  previous_rows <- df[df$time == previous_tw, ]
  
  if (nrow(current_rows) == 0 || nrow(previous_rows) == 0) next
  
  for (j in 1:nrow(current_rows)) {
    c_row <- current_rows[j, ]
    c_lin <- c_row$lineage
    
    matched_parent_id <- NA_integer_
    
    # -------------------------------------------------------------
    # STRATEGIA 1: Standardowe dopasowanie dokładne i skracanie od końca
    # -------------------------------------------------------------
    search_str <- c_lin
    while (nchar(search_str) > 0) {
      matching_parents <- previous_rows[previous_rows$lineage == search_str, ]
      if (nrow(matching_parents) > 0) {
        matched_parent_id <- matching_parents$row_id[1]
        break
      }
      search_str <- substr(search_str, 1, nchar(search_str) - 1)
    }
    
    # -------------------------------------------------------------
    # STRATEGIA 2: Zaawansowany Regex Fallback (Wybór rodzeństwa wg klucza alfabetycznego)
    # -------------------------------------------------------------
    if (is.na(matched_parent_id)) {
      search_str_regex <- c_lin
      
      while (nchar(search_str_regex) > 0) {
        base_str <- substr(search_str_regex, 1, nchar(search_str_regex) - 1)
        last_char <- substr(search_str_regex, nchar(search_str_regex), nchar(search_str_regex))
        
        pattern_any <- paste0("^", base_str, ".$")
        matching_candidates <- previous_rows[grepl(pattern_any, previous_rows$lineage), ]
        
        if (nrow(matching_candidates) > 0) {
          if (grepl("^[A-Za-z]$", last_char)) {
            matching_candidates$cand_last_char <- substr(matching_candidates$lineage, nchar(matching_candidates$lineage), nchar(matching_candidates$lineage))
            lesser_candidates <- matching_candidates[matching_candidates$cand_last_char < last_char, ]
            
            if (nrow(lesser_candidates) > 0) {
              lesser_candidates <- lesser_candidates[order(lesser_candidates$lineage, decreasing = TRUE), ]
              matched_parent_id <- lesser_candidates$row_id[1]
              break
            }
          }
          
          matching_candidates <- matching_candidates[order(matching_candidates$lineage), ]
          matched_parent_id <- matching_candidates$row_id[1]
          break
        }
        search_str_regex <- substr(search_str_regex, 1, nchar(search_str_regex) - 1)
      }
    }
    
    # Przypisujemy ID odnalezionego rodzica do głównej tabeli
    df$parent_id[df$row_id == c_row$row_id] <- matched_parent_id
  }
}

# 2.5. Tabela powiązań rodziców (MAPPING TSV)
parents_mapping <- df %>%
  left_join(df, by = c("parent_id" = "row_id"), suffix = c("_child", "_parent")) %>%
  select(
    time_window        = time_child,
    population         = population_child,
    lineage            = lineage_child,
    parent_time_window = time_parent,
    parent_population  = population_parent,
    parent_lineage     = lineage_parent
  ) %>%
  mutate(
    parent_time_window = ifelse(is.na(parent_time_window), "ROOT", parent_time_window),
    parent_population  = ifelse(is.na(parent_population), "NONE", parent_population),
    parent_lineage     = ifelse(is.na(parent_lineage), "NONE", parent_lineage)
  )

write.table(parents_mapping, output_file_parents, sep = "\t", row.names = FALSE, quote = FALSE)


# 3. Funkcja rekurencyjna zapisująca wcięte drzewa tekstowe (Zadanie stare)
print_trees_recursive <- function(current_id, indent, conn_full, conn_lin) {
  node <- df[df$row_id == current_id, ]
  
  line_full <- paste0(indent, node$time, " ", node$population, ", ", node$lineage)
  writeLines(line_full, conn_full)
  
  line_lin <- paste0(indent, node$time, " ", node$lineage)
  writeLines(line_lin, conn_lin)
  
  children <- df[which(df$parent_id == current_id), ]
  
  if (nrow(children) > 0) {
    children <- children[order(children$lineage), ]
    for (i in 1:nrow(children)) {
      print_trees_recursive(children$row_id[i], paste0(indent, "  "), conn_full, conn_lin)
    }
  }
}

# 3.5. NOWOŚĆ: Funkcja rekurencyjna do generowania pełnych ścieżek od korzenia do liści
find_leaf_paths_recursive <- function(current_id, current_path_vector, conn_paths) {
  node <- df[df$row_id == current_id, ]
  
  # Formatujemy aktualny węzeł jako element krotki: (time_window;population;lineage)
  node_string <- paste0("(", node$time, ";", node$population, ";", node$lineage, ")")
  
  # Dokładamy bieżący element do dotychczas zebranej ścieżki pionowej
  updated_path_vector <- c(current_path_vector, node_string)
  
  # Szukamy dzieci dla tego klastra
  children <- df[which(df$parent_id == current_id), ]
  
  if (nrow(children) == 0) {
    # LIŚĆ: Jeśli klaster nie ma dzieci, ścieżka się kończy. Łączymy elementy separatorem "|" i zapisujemy
    full_path_line <- paste(updated_path_vector, collapse = "|")
    writeLines(full_path_line, conn_paths)
  } else {
    # WĘZEŁ WEWNĘTRZNY: Idziemy głębiej we wszystkie gałęzie potomne (posortowane alfabetycznie)
    children <- children[order(children$lineage), ]
    for (i in 1:nrow(children)) {
      find_leaf_paths_recursive(children$row_id[i], updated_path_vector, conn_paths)
    }
  }
}


# 4. Otwarcie połączeń i jednoczesny zapis struktur danych
file_full  <- file(output_file_full, open = "w")
file_lin   <- file(output_file_lineage, open = "w")
file_paths <- file(output_file_paths, open = "w")

# Wyłonienie klastrów bazowych (korzeni drzewa)
roots <- df[is.na(df$parent_id), ]
roots <- roots[order(roots$time, roots$lineage), ]

# Uruchomienie obu algorytmów rekurencyjnych dla każdego korzenia
for (i in 1:nrow(roots)) {
  # Zapis drzew wciętych
  print_trees_recursive(roots$row_id[i], "", file_full, file_lin)
  
  # Zapis ścieżek od korzenia do liścia (zaczynamy z pustym wektorem ścieżki)
  find_leaf_paths_recursive(roots$row_id[i], character(), file_paths)
}

# Bezpieczne zamknięcie wszystkich strumieni plików
close(file_full)
close(file_lin)
close(file_paths)

# Komunikat końcowy
cat("Sukces! Wszystkie cztery pliki zostały pomyślnie wygenerowane:\n",
    "1. Pełne drzewo tekstowe:      ", output_file_full, "\n",
    "2. Drzewo samych linii:         ", output_file_lineage, "\n",
    "3. Tabela powiązań (TSV):       ", output_file_parents, "\n",
    "4. Ścieżki od korzenia do liścia:", output_file_paths, "\n")
