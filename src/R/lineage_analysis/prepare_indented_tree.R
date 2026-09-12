library(dplyr)
library(stringr)

config_path <- file.path("config", "config.yml")
config <- yaml::yaml.load_file(config_path)

resultsdir <- config$paths$resultsdir
input_file <- file.path(resultsdir, "lineage_analysis", "lineage_populations.tsv")

output_file_full    <- file.path(resultsdir, "lineage_analysis", "cluster_tree_indented.txt")
output_file_lineage <- file.path(resultsdir, "lineage_analysis", "lineage_only_tree_indented.txt")
output_file_parents <- file.path(resultsdir, "lineage_analysis", "lineage_parents_mapping.tsv")
output_file_paths   <- file.path(resultsdir, "lineage_analysis", "root_to_leaf_paths.txt")

df <- read.table(input_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
df$row_id <- 1:nrow(df)
df$parent_id <- NA_integer_

time_windows <- c("00-02", "02-04", "04-06", "06-08", "08-10", 
                  "10-12", "12-14", "14-16", "16-18", "18-20")

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
    
    search_str <- c_lin
    while (nchar(search_str) > 0) {
      matching_parents <- previous_rows[previous_rows$lineage == search_str, ]
      if (nrow(matching_parents) > 0) {
        matched_parent_id <- matching_parents$row_id[1]
        break
      }
      search_str <- substr(search_str, 1, nchar(search_str) - 1)
    }
    
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
    
    df$parent_id[df$row_id == c_row$row_id] <- matched_parent_id
  }
}

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

find_leaf_paths_recursive <- function(current_id, current_path_vector, conn_paths) {
  node <- df[df$row_id == current_id, ]
  
  node_string <- paste0("(", node$time, ";", node$population, ";", node$lineage, ")")
  
  updated_path_vector <- c(current_path_vector, node_string)
  

  children <- df[which(df$parent_id == current_id), ]
  
  if (nrow(children) == 0) {
    full_path_line <- paste(updated_path_vector, collapse = "|")
    writeLines(full_path_line, conn_paths)
  } else {
    children <- children[order(children$lineage), ]
    for (i in 1:nrow(children)) {
      find_leaf_paths_recursive(children$row_id[i], updated_path_vector, conn_paths)
    }
  }
}


file_full  <- file(output_file_full, open = "w")
file_lin   <- file(output_file_lineage, open = "w")
file_paths <- file(output_file_paths, open = "w")

roots <- df[is.na(df$parent_id), ]
roots <- roots[order(roots$time, roots$lineage), ]

for (i in 1:nrow(roots)) {
  print_trees_recursive(roots$row_id[i], "", file_full, file_lin)
  
  find_leaf_paths_recursive(roots$row_id[i], character(), file_paths)
}

close(file_full)
close(file_lin)
close(file_paths)

cat("Success! Files:\n",
    "1. Tree of populations:      ", output_file_full, "\n",
    "2. Tree of lineages:         ", output_file_lineage, "\n",
    "3. Parents table (TSV):       ", output_file_parents, "\n",
    "4. Paths from root to leaves:", output_file_paths, "\n")
