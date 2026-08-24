library(data.table)

config_path <- file.path("config", "config.yml")
config <- yaml::yaml.load_file(config_path)

datadir <- config$paths$datadir
resultsdir <- config$paths$resultsdir

atac_meta_path <- file.path(datadir, "calderon_data", "atac_meta.rds")
nn_dir <- file.path(datadir, "calderon_data", "new_time", "NN")

output_file_gz <- file.path(resultsdir, "cluster_cardinality.tsv.gz")

atac_meta <- readRDS(atac_meta_path)
class(atac_meta)
# [1] "data.frame"
head(atac_meta)
#                                                                              cell
# CTATGGTTCGTTCCATTCTTTCTCATTGCCCTCTGCGATC CTATGGTTCGTTCCATTCTTTCTCATTGCCCTCTGCGATC
# AGATAATTCCTTCCATTCTTAGTCGCGTCGAAGTTCGCTG AGATAATTCCTTCCATTCTTAGTCGCGTCGAAGTTCGCTG
# TCCTCTTAACCTAGCTTCTTGTAGTAGTCCTGATTCTCGT TCCTCTTAACCTAGCTTCTTGTAGTAGTCCTGATTCTCGT
# ACGGCAAGCAATATCTTCCGGAGCTCAGCCCCTACTCAAC ACGGCAAGCAATATCTTCCGGAGCTCAGCCCCTACTCAAC
# AACTTACGCTCATTCTGATGCCAATTCCATGTCGTTCGCC AACTTACGCTCATTCTGATGCCAATTCCATGTCGTTCGCC
# CGCGTTACCTTCGGAGTATGGTTCGAGGCGTTAATTCGTA CGCGTTACCTTCGGAGTATGGTTCGAGGCGTTAATTCGTA
#                                          doublet_score  exp  time
# CTATGGTTCGTTCCATTCTTTCTCATTGCCCTCTGCGATC     0.8027211 exp1 03-07
# AGATAATTCCTTCCATTCTTAGTCGCGTCGAAGTTCGCTG     0.7006803 exp1 03-07
# TCCTCTTAACCTAGCTTCTTGTAGTAGTCCTGATTCTCGT     0.5051020 exp1 03-07
# ACGGCAAGCAATATCTTCCGGAGCTCAGCCCCTACTCAAC     0.6530612 exp1 03-07
# AACTTACGCTCATTCTGATGCCAATTCCATGTCGTTCGCC     0.6717687 exp1 03-07
# CGCGTTACCTTCGGAGTATGGTTCGAGGCGTTAATTCGTA     0.7789116 exp1 03-07
#                                                    sample seurat_clusters
# CTATGGTTCGTTCCATTCTTTCTCATTGCCCTCTGCGATC exp1_hrs03-07_b1               7
# AGATAATTCCTTCCATTCTTAGTCGCGTCGAAGTTCGCTG exp1_hrs03-07_b1               6
# TCCTCTTAACCTAGCTTCTTGTAGTAGTCCTGATTCTCGT exp1_hrs03-07_b1               1
# ACGGCAAGCAATATCTTCCGGAGCTCAGCCCCTACTCAAC exp1_hrs03-07_b1               1
# AACTTACGCTCATTCTGATGCCAATTCCATGTCGTTCGCC exp1_hrs03-07_b1               1
# CGCGTTACCTTCGGAGTATGGTTCGAGGCGTTAATTCGTA exp1_hrs03-07_b1               1
#                                          lasso_age NNv1_age lasso_time.new
# CTATGGTTCGTTCCATTCTTTCTCATTGCCCTCTGCGATC -3.110988 1.690188          00-02
# AGATAATTCCTTCCATTCTTAGTCGCGTCGAAGTTCGCTG -1.039036 1.330054          00-02
# TCCTCTTAACCTAGCTTCTTGTAGTAGTCCTGATTCTCGT -3.924070 1.984080          00-02
# ACGGCAAGCAATATCTTCCGGAGCTCAGCCCCTACTCAAC  2.537512 1.861791          02-04
# AACTTACGCTCATTCTGATGCCAATTCCATGTCGTTCGCC -1.401258 1.801240          00-02
# CGCGTTACCTTCGGAGTATGGTTCGAGGCGTTAATTCGTA  2.353294 1.718721          02-04
#                                          NNv1_time.new seurat_clusters.predtime
# CTATGGTTCGTTCCATTCTTTCTCATTGCCCTCTGCGATC         00-02                        3
# AGATAATTCCTTCCATTCTTAGTCGCGTCGAAGTTCGCTG         00-02                        1
# TCCTCTTAACCTAGCTTCTTGTAGTAGTCCTGATTCTCGT         00-02                        1
# ACGGCAAGCAATATCTTCCGGAGCTCAGCCCCTACTCAAC         00-02                        2
# AACTTACGCTCATTCTGATGCCAATTCCATGTCGTTCGCC         00-02                        2
# CGCGTTACCTTCGGAGTATGGTTCGAGGCGTTAATTCGTA         00-02                        2
#                                          refined_annotation
# CTATGGTTCGTTCCATTCTTTCTCATTGCCCTCTGCGATC         Blastoderm
# AGATAATTCCTTCCATTCTTAGTCGCGTCGAAGTTCGCTG         Blastoderm
# TCCTCTTAACCTAGCTTCTTGTAGTAGTCCTGATTCTCGT         Blastoderm
# ACGGCAAGCAATATCTTCCGGAGCTCAGCCCCTACTCAAC         Blastoderm
# AACTTACGCTCATTCTGATGCCAATTCCATGTCGTTCGCC         Blastoderm
# CGCGTTACCTTCGGAGTATGGTTCGAGGCGTTAATTCGTA         Blastoderm

print(unique(atac_meta$NNv1_time.new))
#  [1] 00-02 02-04 04-06 06-08 08-10 10-12 12-14 14-16 16-18 18-20
# Levels: 00-02 02-04 04-06 06-08 08-10 10-12 12-14 14-16 16-18 18-20

print(unique(atac_meta$seurat_clusters.predtime))
#  [1] 3  1  2  4  0  5  6  7  8  11 18 14 10 9  17 12 15 16 13 19 21 20 25 24 23
# [26] 22 26
# 27 Levels: 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 ... 26

print(unique(atac_meta$refined_annotation))
#  [1] Blastoderm               Unknown                  Ectoderm anlage         
#  [4] Endoderm anlage          Mesoderm anlage          Germ cell               
#  [7] Muscle prim.             Hindgut prim.            Yolk                    
# [10] Midgut prim.             Head ectoderm            Plasmatocytes           
# [13] Amnioserosa anlage       Ventral nerve cord prim. Mesectoderm anlage      
# [16] Foregut prim.            Epidermis prim.          Head ectoderm prim.     
# [19] Tracheal system prim.    PNS & sense              Brain prim.             
# [22] Amnioserosa              Ventral midline          Fat body                
# [25] Glia                     Ventral nerve cord       Somatic muscle          
# [28] Brain                    Epidermis                Visceral muscle         
# [31] Midgut                   Tracheal system          Hindgut                 
# [34] Pharnyx                  Malpighian tubule        Salivary gland          
# [37] Neural                   Proventriculus          
# 38 Levels: Amnioserosa Amnioserosa anlage Blastoderm Brain ... Yolk

# from atac_meta.readme.txt
# column_name,notes
# NNv1_time.new,NN model predicted 2 hr time window
# seurat_clusters.predtime,Cluster id based on clustering in each NN model predicted 2hr time window
# refined_annotation,Annotation based on each cluster id based on clustering in each NN model predicted 2hr time window

cluster_cardinality <- table(atac_meta$seurat_clusters.predtime, atac_meta$NNv1_time.new)
print(cluster_cardinality)
#      00-02 02-04 04-06 06-08 08-10 10-12 12-14 14-16 16-18 18-20
#   0  11776 22071 15229 39573 29461 18471 23524 12976  4865   280
#   1   8128 16391 15103 32325 19587 11230 15156 11551  4217   234
#   2   6967  8777 14968 29463 19496  9804 13584 10280  3976   181
#   3   2058  7702 12728 15815 13648  8164 13009  9795  3846   164
#   4   1255  5124 11065 15245 11617  7894  9095  8523  3503   132
#   5    835  3322 10334 14440  5082  7673  7376  6760  3453   114
#   6      0  2015  6994 13456  4867  6793  7182  6427  3363    90
#   7      0   989  6821  5880  4519  6488  7064  5385  2687    87
#   8      0   172  6637  4645  4363  6066  7016  4922  2632    75
#   9      0     0  5226  4262  3917  4984  6666  4905  2363    67
#   10     0     0  4458  4137  3623  4842  6540  4876  2010     0
#   11     0     0  4618  3597  2719  4529  5398  4669  1917     0
#   12     0     0  3816  3111  2593  3814  3794  3934  1518     0
#   13     0     0  2301  3009  2518  3259  3246  3761  1315     0
#   14     0     0  2289  2026  2420  3220  3176  2744  1273     0
#   15     0     0  1648  1031  1763  2167  2195  2186  1262     0
#   16     0     0  1023   765   629  1345  1944  2059  1165     0
#   17     0     0  1019     0   617  1315  1815  2017  1157     0
#   18     0     0   792     0   560  1000  1560  1978  1004     0
#   19     0     0     0     0     0   694  1149  1874   903     0
#   20     0     0     0     0     0   635   942  1407   873     0
#   21     0     0     0     0     0   544   676   714   745     0
#   22     0     0     0     0     0     0   667   247   400     0
#   23     0     0     0     0     0     0   136   217   202     0
#   24     0     0     0     0     0     0    96   180   189     0
#   25     0     0     0     0     0     0    95   126   144     0
#   26     0     0     0     0     0     0     0    75     0     0

#######################################################################################

dt_output <- as.data.table(as.data.frame.matrix(cluster_cardinality), keep.rownames = "seurat_clusters_predtime")
fwrite(dt_output, file = output_file_gz, sep = "\t", compress = "gzip")

# The resulting file can be later read using:
#card <- fread(output_file_gz)
#print(card)

#   seurat_clusters_predtime 00-02 02-04 04-06 06-08 08-10 10-12 12-14 14-16
#                       <int> <int> <int> <int> <int> <int> <int> <int> <int>
# 1:                        0 11776 22071 15229 39573 29461 18471 23524 12976
# 2:                        1  8128 16391 15103 32325 19587 11230 15156 11551
# 3:                        2  6967  8777 14968 29463 19496  9804 13584 10280
# 4:                        3  2058  7702 12728 15815 13648  8164 13009  9795
# 5:                        4  1255  5124 11065 15245 11617  7894  9095  8523
# 6:                        5   835  3322 10334 14440  5082  7673  7376  6760
# 7:                        6     0  2015  6994 13456  4867  6793  7182  6427
# 8:                        7     0   989  6821  5880  4519  6488  7064  5385
# 9:                        8     0   172  6637  4645  4363  6066  7016  4922
#10:                        9     0     0  5226  4262  3917  4984  6666  4905
#11:                       10     0     0  4458  4137  3623  4842  6540  4876
#12:                       11     0     0  4618  3597  2719  4529  5398  4669
#13:                       12     0     0  3816  3111  2593  3814  3794  3934
#14:                       13     0     0  2301  3009  2518  3259  3246  3761
#15:                       14     0     0  2289  2026  2420  3220  3176  2744
#16:                       15     0     0  1648  1031  1763  2167  2195  2186
#17:                       16     0     0  1023   765   629  1345  1944  2059
#18:                       17     0     0  1019     0   617  1315  1815  2017
#19:                       18     0     0   792     0   560  1000  1560  1978
#20:                       19     0     0     0     0     0   694  1149  1874
#21:                       20     0     0     0     0     0   635   942  1407
#22:                       21     0     0     0     0     0   544   676   714
#23:                       22     0     0     0     0     0     0   667   247
#24:                       23     0     0     0     0     0     0   136   217
#25:                       24     0     0     0     0     0     0    96   180
#26:                       25     0     0     0     0     0     0    95   126
#27:                       26     0     0     0     0     0     0     0    75
#    seurat_clusters_predtime 00-02 02-04 04-06 06-08 08-10 10-12 12-14 14-16
#    16-18 18-20
#    <int> <int>
# 1:  4865   280
# 2:  4217   234
# 3:  3976   181
# 4:  3846   164
# 5:  3503   132
# 6:  3453   114
# 7:  3363    90
# 8:  2687    87
# 9:  2632    75
#10:  2363    67
#11:  2010     0
#12:  1917     0
#13:  1518     0
#14:  1315     0
#15:  1273     0
#16:  1262     0
#17:  1165     0
#18:  1157     0
#19:  1004     0
#20:   903     0
#21:   873     0
#22:   745     0
#23:   400     0
#24:   202     0
#25:   189     0
#26:   144     0
#27:     0     0
#    16-18 18-20

