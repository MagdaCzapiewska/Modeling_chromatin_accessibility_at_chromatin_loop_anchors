library(dplyr)

config_path <- file.path("config", "config.yml")
config <- yaml::yaml.load_file(config_path)

datadir <- config$paths$datadir
resultsdir <- config$paths$resultsdir

input_file <- file.path(resultsdir, "lineage_analysis", "lineage_populations.tsv")
output_file <- file.path(resultsdir, "lineage_analysis", "time_window_clusters.txt")


df <- read.table(input_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)

time_windows <- c("00-02", "02-04", "04-06", "06-08", "08-10", 
                  "10-12", "12-14", "14-16", "16-18", "18-20")

file_conn <- file(output_file, open = "w")

for (tw in time_windows) {

  writeLines(tw, file_conn)
  
  df_subset <- df %>% filter(time == tw)
  

  if (nrow(df_subset) > 0) {
    lines_to_write <- paste(df_subset$population, df_subset$lineage, sep = ", ")
    
    writeLines(lines_to_write, file_conn)
  }
}

close(file_conn)

print(paste("File saved as:", output_file))
