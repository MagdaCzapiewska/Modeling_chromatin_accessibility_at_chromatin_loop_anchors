library(data.table)

config_path <- file.path("config", "config.yml")
config <- yaml::yaml.load_file(config_path)

datadir <- config$paths$datadir
resultsdir <- config$paths$resultsdir

atac_meta_path <- file.path(datadir, "calderon_data", "atac_meta.rds")
output_map_file <- file.path(resultsdir, "cluster_pop_mapping.tsv.gz")

atac_meta <- readRDS(atac_meta_path)
dt <- as.data.table(atac_meta)

pop_map <- unique(dt[, .(
  time_window = as.character(NNv1_time.new),
  cluster_id = as.character(seurat_clusters.predtime),
  refined_annotation
)])

pop_map[, clean_annot := gsub("\\.", "", refined_annotation)]
pop_map[, clean_annot := gsub(" ", "_", clean_annot)]
pop_map[, population := paste0(cluster_id, "_", clean_annot)]

pop_map[, clean_annot := NULL]

setkey(pop_map, time_window, cluster_id)
fwrite(pop_map, file = output_map_file, sep = "\t", compress = "gzip")

message("Mapping saved to: ", output_map_file)
