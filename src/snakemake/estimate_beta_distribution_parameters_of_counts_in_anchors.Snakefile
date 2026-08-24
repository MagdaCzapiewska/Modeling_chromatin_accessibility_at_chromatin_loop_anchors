import os
import pandas as pd

configfile: "config/config.yml"

time_windows = ["00-02", "02-04", "04-06", "06-08", "08-10", "10-12", "12-14", "14-16", "16-18", "18-20"]

RESULTSDIR = config["paths"]["resultsdir"]
DATADIR = config["paths"]["datadir"]

COUNTS_DIR = os.path.join(RESULTSDIR, "counts", "counts_in_anchors")
OUTDIR = os.path.join(RESULTSDIR, "parameters", "beta_distribution_of_counts_in_anchors")
RSCRIPT = os.path.join(config['paths']['Rsrcdir'], "estimate_beta_distribution_parameters_of_counts_in_anchors.R")

loops_file = os.path.join(DATADIR, "long_and_short_range_loops_D_mel.tsv")
loops_table = pd.read_csv(loops_file, sep="\t")
loop_ids = loops_table["loop_id"].dropna().unique().tolist()

##############################################
# Rule all = aim of the pipeline
##############################################
rule all:
    input:
        expand(os.path.join(OUTDIR, "beta_parameters_{tw}_{loop_id}_all_and_by_population.tsv.gz"), tw=time_windows, loop_id=loop_ids)


##############################################
# One time window = one job, one loop
##############################################
rule calculate_reads:
    input:
        script = RSCRIPT,
        reads = os.path.join(COUNTS_DIR, "reads_{tw}_{loop_id}.tsv.gz")
    output:
        tsv = os.path.join(OUTDIR, "beta_parameters_{tw}_{loop_id}_all_and_by_population.tsv.gz")
    log:
        os.path.join(OUTDIR, "logs", "{tw}_{loop_id}.log")
    shell:
        """
        Rscript {input.script} "{wildcards.tw}" "{wildcards.loop_id}" > {log} 2>&1
        """
