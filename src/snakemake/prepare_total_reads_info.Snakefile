import os

configfile: "config/config.yml"

time_windows = ["00-02", "02-04", "04-06", "06-08", "08-10", "10-12", "12-14", "14-16", "16-18", "18-20"]

RESULTSDIR = config["paths"]["resultsdir"]
DATADIR = config["paths"]["datadir"]

OUTDIR = os.path.join(RESULTSDIR, "counts", "total_reads")
RSCRIPT = os.path.join(config['paths']['Rsrcdir'], "prepare_total_reads_info.R")

##############################################
# Rule all = aim of the pipeline
##############################################
rule all:
    input:
        expand(os.path.join(OUTDIR, "total_reads_{tw}.tsv.gz"), tw=time_windows)


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
        tsv_gz = os.path.join(OUTDIR, "total_reads_{tw}.tsv.gz")
    shell:
        """
        Rscript {input.script} "{wildcards.tw}"
        """
