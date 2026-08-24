import os
import pandas as pd

configfile: "config/config.yml"

time_windows = ["00-02", "02-04", "04-06", "06-08", "08-10", "10-12", "12-14", "14-16", "16-18", "18-20"]

RESULTSDIR = config["paths"]["resultsdir"]
DATADIR = config["paths"]["datadir"]

OUTDIR = os.path.join(RESULTSDIR, "counts", "counts_in_anchors")
RSCRIPT = os.path.join(config['paths']['Rsrcdir'], "prepare_counts_in_anchors.R")

loops_file = os.path.join(DATADIR, "long_and_short_range_loops_D_mel.tsv")
loops_table = pd.read_csv(loops_file, sep="\t")
loop_ids = loops_table["loop_id"].dropna().unique().tolist()

##############################################
# Rule all = aim of the pipeline
##############################################
rule all:
    input:
        expand(os.path.join(OUTDIR, "reads_{tw}_{loop}.tsv.gz"), tw=time_windows, loop=loop_ids)


##############################################
# One time window = one job

# Before running shell section, Snakemake automatically creates
# directories for all files listed in output section.
##############################################
rule calculate_reads:
    input:
        script = RSCRIPT,
        meta = os.path.join(DATADIR, "calderon_data", "atac_meta.rds"),
        matrix = os.path.join(DATADIR, "calderon_data", "new_time", "NN", "GSE190130_{tw}_new_timeNN.peak_matrix.mtx.gz"),
        peaks = os.path.join(DATADIR, "calderon_data", "new_time", "NN", "GSE190130_{tw}_new_timeNN.peak_matrix.rows.txt.gz"),
        barcodes = os.path.join(DATADIR, "calderon_data", "new_time", "NN", "GSE190130_{tw}_new_timeNN.peak_matrix.columns.txt.gz")
    output:
        # We inform Snakemake, that this single job creates files for all loops in this time window
        tsv_gz = expand(os.path.join(OUTDIR, "reads_{{tw}}_{loop}.tsv.gz"), loop=loop_ids)
    shell:
        """
        Rscript {input.script} "{wildcards.tw}"
        """
