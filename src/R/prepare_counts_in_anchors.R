library(GenomicRanges)
library(Seurat)
library(hdf5r)
library(Matrix)
library(data.table)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
  stop("Missing required time window argument (tw). Please provide e.g. '00-02'.")
}

tw <- args[1]

config_path <- file.path("config", "config.yml")
config <- yaml::yaml.load_file(config_path)

datadir <- config$paths$datadir
resultsdir <- config$paths$resultsdir
nn_dir <- file.path(datadir, "calderon_data", "new_time", "NN")
total_reads_dir <- file.path(resultsdir, "counts", "total_reads")

atac_meta_path <- file.path(datadir, "calderon_data", "atac_meta.rds")
loops_path <- file.path(datadir, "long_and_short_range_loops_D_mel.tsv")

out_dir <- file.path(resultsdir, "counts", "counts_in_anchors")
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

################################################################################

atac_meta <- readRDS(atac_meta_path)
loops <- read.delim(loops_path, header = TRUE, sep = "\t")

stopifnot(length(unique(loops$loop_id)) == nrow(loops))
dim(loops)
# [1] 417  20
loops_nona <- loops[!is.na(loops$x1) & !is.na(loops$x2) & 
                               !is.na(loops$y1) & !is.na(loops$y2), ]
dim(loops_nona)
# [1] 417  20

loops_subset <- loops_nona

left_anchors <- GRanges(seqnames = loops_subset$chr1,
                        ranges = IRanges(start = loops_subset$x1,
                                         end = loops_subset$x2))

right_anchors <- GRanges(seqnames = loops_subset$chr2,
                         ranges = IRanges(start = loops_subset$y1,
                                          end = loops_subset$y2))

loop_ids <- loops_subset$loop_id

loop_anchors <- list(
  left = left_anchors,
  right = right_anchors,
  loop_ids = loop_ids
)

################################################################################

get_anchor_peak_hits <- function(peaks_gr, loop_anchors) {
  return(list(
    hits1 = findOverlaps(loop_anchors$left, peaks_gr),
    hits2 = findOverlaps(loop_anchors$right, peaks_gr)
  ))
}

peaks_to_granges <- function(peak_names) {
  parsed <- do.call(rbind, strsplit(peak_names, "_"))
  chr_raw <- parsed[, 1]
  if (any(!grepl("^chr", chr_raw))) stop("All chromosome names must start with 'chr'")
  chr <- sub("^chr", "", chr_raw)
  bed_start <- as.integer(parsed[, 2])
  bed_end <- as.integer(parsed[, 3])
  
  return(GRanges(
    seqnames = chr,
    ranges = IRanges(start = bed_start + 1, end = bed_end),
    names = peak_names
  ))
}

##################################################################################

total_reads_file <- file.path(total_reads_dir, paste0("total_reads_", tw, ".tsv.gz"))
  
if (!file.exists(total_reads_file)) {
  stop(paste("Required total reads file does not exist for time window:", tw))
}
  
cat("\n=== Processing time window:", tw, "===\n")
total_reads_df <- fread(total_reads_file)
  
prefix <- paste0("GSE190130_", tw, "_new_timeNN.peak_matrix")
columns_file <- file.path(nn_dir, paste0(prefix, ".columns.txt.gz"))
rows_file <- file.path(nn_dir, paste0(prefix, ".rows.txt.gz"))
matrix_file <- file.path(nn_dir, paste0(prefix, ".mtx.gz"))
  
barcodes <- readLines(gzfile(columns_file))
peaks <- readLines(gzfile(rows_file))
mat <- readMM(gzfile(matrix_file))

# Conversion from dgTMatrix to dgCMatrix for safe indexing and quick colSums
mat <- as(mat, "CsparseMatrix")

stopifnot(setequal(barcodes, total_reads_df$barcode))
  
original_order <- total_reads_df$barcode
reordering_index <- match(barcodes, original_order)
total_reads_df <- total_reads_df[reordering_index, ]

if (!all(reordering_index == seq_along(barcodes))) {
  cat("There was a change in row order in total_reads_df!\n")
}
stopifnot(all(reordering_index == seq_along(barcodes)))

cat("Time window, ", tw, ", Start peaks_to_granges\n")
peaks_gr <- peaks_to_granges(peaks)
cat("Time window, ", tw, ", Stop peaks_to_granges\n")
    
anchors <- loop_anchors
num_loops <- length(anchors$left)
cat("  >> Number of loops:", num_loops, "\n")

anchor_hits <- get_anchor_peak_hits(peaks_gr, loop_anchors)
cat("Time window, ", tw, ", After get_anchor_peak_hits\n")

loop_to_peaks_left <- split(
  subjectHits(anchor_hits$hits1),
  factor(queryHits(anchor_hits$hits1), levels = seq_len(num_loops)),
  drop = FALSE
)

loop_to_peaks_right <- split(
  subjectHits(anchor_hits$hits2),
  factor(queryHits(anchor_hits$hits2), levels = seq_len(num_loops)),
  drop = FALSE
)
cat("After split\n")

for (loop_idx in seq_len(num_loops)) {
  cat("Time window ", tw, ", ", loop_idx, "/", num_loops, "loops processed\n")
  
  left_peak_idx <- loop_to_peaks_left[[loop_idx]]
  right_peak_idx <- loop_to_peaks_right[[loop_idx]]
    
  reads_A1 <- if (length(left_peak_idx) > 0) {
    as.numeric(colSums(mat[left_peak_idx, , drop = FALSE]))
  } else {
    rep(0, ncol(mat))
  }
    
  reads_A2 <- if (length(right_peak_idx) > 0) {
    as.numeric(colSums(mat[right_peak_idx, , drop = FALSE]))
  } else {
    rep(0, ncol(mat))
  }
    
  loop_id <- anchors$loop_ids[loop_idx]
  
  out_loop_dt <- data.table(
    barcode = total_reads_df$barcode,
    population = total_reads_df$population,
    total_reads = total_reads_df$total_reads,
    x_A1 = reads_A1,
    x_A2 = reads_A2
  )
  
  output_file <- file.path(out_dir, paste0("reads_", tw, "_", loop_id, ".tsv.gz"))
  fwrite(out_loop_dt, file = output_file, sep = "\t", quote = FALSE, compress = "gzip")
}
  
cat("Ended time window:", tw, "\n")
