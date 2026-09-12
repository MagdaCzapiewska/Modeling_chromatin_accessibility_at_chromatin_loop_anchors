library(readxl)
library(dplyr)
library(stringr)

config_path <- file.path("config", "config.yml")
config <- yaml::yaml.load_file(config_path)

datadir <- config$paths$datadir
resultsdir <- config$paths$resultsdir

input_file <- file.path(datadir, "calderon_data", "science.abn5800_table_s2.xlsx")


outdir <- file.path(resultsdir, "lineage_analysis")

if (!dir.exists(outdir)) {
  dir.create(outdir, recursive = TRUE)
}

df <- read_excel(input_file, skip = 1)

df_processed <- df %>%
  select(time, cluster, lineage, refined_annotation) %>%
  mutate(
    clean_annotation = str_replace_all(refined_annotation, "\\.", ""),
    clean_annotation = str_replace_all(clean_annotation, " ", "_"),
    
    population = paste(cluster, clean_annotation, sep = "_")
  ) %>%
  select(time, cluster, lineage, refined_annotation, population)

output_file <- file.path(outdir, "lineage_populations.tsv")

write.table(df_processed, file = output_file, sep = "\t", quote = FALSE, row.names = FALSE)

print(paste("File saved as:", output_file))
