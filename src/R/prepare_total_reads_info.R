library(data.table)
library(Matrix)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
  stop("Missing required time window argument (tw). Please provide e.g. '00-02'.")
}

tw <- args[1]

config_path <- file.path("config", "config.yml")
config <- yaml::yaml.load_file(config_path)

datadir <- config$paths$datadir
resultsdir <- config$paths$resultsdir

atac_meta_path <- file.path(datadir, "calderon_data", "atac_meta.rds")
nn_dir <- file.path(datadir, "calderon_data", "new_time", "NN")

out_dir <- file.path(resultsdir, "counts", "total_reads")
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}
out_file <- file.path(out_dir, paste0("total_reads_", tw, ".tsv.gz"))

atac_meta <- readRDS(atac_meta_path)

###################################### Inspection of the data ###########################################

# prefix <- file.path(nn_dir, "GSE190130_00-02_new_timeNN.peak_matrix")

# columns_file <- paste0(prefix, ".columns.txt.gz")
# rows_file <- paste0(prefix, ".rows.txt.gz")
# matrix_file <- paste0(prefix, ".mtx.gz")

# cell_barcodes <- readLines(gzfile(columns_file))
# peak_names <- readLines(gzfile(rows_file))
# peak_matrix <- readMM(gzfile(matrix_file))

# length(cell_barcodes)
# [1] 31019

# length(peak_names)
# [1] 110185

# head(cell_barcodes)
# [1] "AACTATTGGCCGCGCAACCGGATAGCCGATTCTAGCCTCC"
# [2] "AACTTACGCTGTACCTGAATTAAGTAAGTAACGTTAAGAC"
# [3] "AAGCCGGTAAACGAGGAGCCCTGATCAGAGGTCGCGAATA"
# [4] "AAGCCGGTAAGGCGGATACCGCGCCTAGTTACGTTAAGAC"
# [5] "AAGGACCTTAATATCTTCCGCTGATCAGAGATCTCCAGCT"
# [6] "AAGTCTTCTGCTGCATTACGCCAATTCCATGCATACCAAT"

# head(peak_names)
# [1] "chr2L_5534_5763"   "chr2L_5770_6000"   "chr2L_6586_6876"  
# [4] "chr2L_7483_9420"   "chr2L_12086_12236" "chr2L_12396_13010"

# dim(peak_matrix)
# [1] 110185  31019

# colnames(peak_matrix) <- cell_barcodes
# rownames(peak_matrix) <- peak_names

# dim(peak_matrix)  # number of peaks x number of cells
# [1] 110185  31019

# str(peak_matrix)
# Formal class 'dgTMatrix' [package "Matrix"] with 6 slots
#   ..@ i       : int [1:76433469] 0 0 0 0 0 0 0 0 0 0 ...
#   ..@ j       : int [1:76433469] 14 58 179 203 210 253 255 274 279 290 ...
#   ..@ Dim     : int [1:2] 110185 31019
#   ..@ Dimnames:List of 2
#   .. ..$ : chr [1:110185] "chr2L_5534_5763" "chr2L_5770_6000" "chr2L_6586_6876" "chr2L_7483_9420" ...
#   .. ..$ : chr [1:31019] "AACTATTGGCCGCGCAACCGGATAGCCGATTCTAGCCTCC" "AACTTACGCTGTACCTGAATTAAGTAAGTAACGTTAAGAC" "AAGCCGGTAAACGAGGAGCCCTGATCAGAGGTCGCGAATA" "AAGCCGGTAAGGCGGATACCGCGCCTAGTTACGTTAAGAC" ...
#   ..@ x       : num [1:76433469] 2 1 1 2 2 1 2 1 2 2 ...
#   ..@ factors : list()

# print(peak_matrix[1:10, 1:10])
# 10 x 10 sparse Matrix of class "dgTMatrix"
#   [[ suppressing 10 column names ‘AACTATTGGCCGCGCAACCGGATAGCCGATTCTAGCCTCC’, ‘AACTTACGCTGTACCTGAATTAAGTAAGTAACGTTAAGAC’, ‘AAGCCGGTAAACGAGGAGCCCTGATCAGAGGTCGCGAATA’ ... ]]
                                     
# chr2L_5534_5763   . . . . . . . . . .
# chr2L_5770_6000   . . . . . . . . . .
# chr2L_6586_6876   . . . . . . . . . .
# chr2L_7483_9420   . . . . . . . . . .
# chr2L_12086_12236 1 . . . . . . . . .
# chr2L_12396_13010 . . . . . . . . . .
# chr2L_13533_14305 . . . . . . . . . .
# chr2L_14323_14895 . . . 1 . . 1 . . .
# chr2L_14966_15273 . . . . . . 1 . . .
# chr2L_15382_15667 . . . . . . . . . .

#########################################################################################################################

prefix <- paste0("GSE190130_", tw, "_new_timeNN.peak_matrix")
columns_file <- file.path(nn_dir, paste0(prefix, ".columns.txt.gz"))
rows_file <- file.path(nn_dir, paste0(prefix, ".rows.txt.gz"))
matrix_file <- file.path(nn_dir, paste0(prefix, ".mtx.gz"))
  
barcodes <- readLines(gzfile(columns_file))
peaks <- readLines(gzfile(rows_file))
mat <- readMM(gzfile(matrix_file))

if (nrow(mat) != length(peaks) || ncol(mat) != length(barcodes)) {
  stop("Matrix dimensions do not match the number of rows or columns.")
}
  
meta_sub <- atac_meta[atac_meta$NNv1_time.new == tw, ]
intersecting_barcodes <- meta_sub$cell[meta_sub$cell %in% barcodes]

if (length(intersecting_barcodes) != length(barcodes)) {
  stop(sprintf("Barcode count mismatch: %d in file vs %d in metadata for window %s.",
               length(barcodes), length(intersecting_barcodes), tw))
}

if (!all(sort(intersecting_barcodes) == sort(barcodes))) {
  stop(sprintf("Barcode set mismatch in time window %s!", tw))
}

meta_sub <- meta_sub[meta_sub$cell %in% barcodes, ]
meta_sub <- meta_sub[match(barcodes, meta_sub$cell), ]

cluster <- meta_sub$seurat_clusters.predtime

annotation <- meta_sub$refined_annotation
annotation <- gsub(" ", "_", annotation)
annotation <- gsub("\\.", "", annotation)

population <- paste(cluster, annotation, sep = "_")

total_reads <- colSums(mat)
number_of_peaks <- colSums(mat > 0)

result_df <- data.frame(
  barcode = barcodes,
  cluster = cluster,
  annotation = annotation,
  population = population,
  total_reads = total_reads,
  number_of_peaks = number_of_peaks,
  stringsAsFactors = FALSE
)

fwrite(result_df, file = out_file, sep = "\t", quote = FALSE, compress = "gzip")

message("Finished counting total reads for time window: ", tw)
