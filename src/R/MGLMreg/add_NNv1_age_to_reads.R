library(data.table)

config_path <- file.path("config", "config.yml")
config <- yaml::yaml.load_file(config_path)

datadir <- config$paths$datadir
resultsdir <- config$paths$resultsdir

input_dir <- file.path(resultsdir, "counts", "counts_in_anchors")
output_dir <- file.path(resultsdir, "counts", "counts_in_anchors")

loops_file <- file.path(datadir, "long_and_short_range_loops_D_mel.tsv")
metadata_file <- file.path(datadir, "calderon_data", "atac_meta.rds")

time_windows <- c("00-02", "02-04", "04-06", "06-08", "08-10", "10-12", "12-14", "12-14", "14-16", "16-18", "18-20")

loops_table <- fread(loops_file)
loops_of_interest <- loops_table$loop_id

metadata <- readRDS(metadata_file)
metadata <- as.data.table(metadata)

setkey(metadata, NNv1_time.new, cell)

for (tw in time_windows) {
    cat("=== Time window:", tw, "===\n")
    
    meta_tw <- metadata[.(tw), .(cell, NNv1_age)]
    
    if (nrow(meta_tw) == 0) {
        cat("  No cells for time window ", tw, "- skipping.\n")
        next
    }
    
    for (loop_id in loops_of_interest) {
        input_file  <- file.path(input_dir, paste0("reads_", tw, "_", loop_id, ".tsv.gz"))
        output_file <- file.path(output_dir, paste0("reads_", tw, "_", loop_id, ".tsv.gz"))

        if (!file.exists(input_file)) {
            next
        }
        
        reads_dt <- fread(input_file)
        
        if (nrow(reads_dt) == 0) {
            next
        }
        
        reads_dt[meta_tw, NNv1_age := i.NNv1_age, on = .(barcode = cell)]

        fwrite(reads_dt, output_file, sep = "\t", compress = "gzip")
    }
}

cat("=== Finished ===\n")
