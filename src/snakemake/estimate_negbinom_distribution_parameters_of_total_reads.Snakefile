import os

configfile: "config/config.yml"

time_windows = ["00-02", "02-04", "04-06", "06-08", "08-10", "10-12", "12-14", "14-16", "16-18", "18-20"]

RESULTSDIR = config["paths"]["resultsdir"]
DATADIR = config["paths"]["datadir"]

TOTAL_READS_DIR = os.path.join(RESULTSDIR, "counts", "total_reads")
OUTDIR = os.path.join(RESULTSDIR, "parameters", "negbinom_distribution_of_total_reads")
RSCRIPT = os.path.join(config['paths']['Rsrcdir'], "estimate_negbinom_distribution_parameters_of_total_reads.R")

##############################################
# Rule all = aim of the project
##############################################
rule all:
    input:
        expand(os.path.join(OUTDIR, "negbinom_parameters_{tw}_all_and_by_population.tsv.gz"), tw=time_windows),
        expand(os.path.join(OUTDIR, "negbinom_models_{tw}_all_and_by_population.rds"), tw=time_windows)


##############################################
# one job = one time window
##############################################
rule estimate_negbinom_parameters:
    input:
        script = RSCRIPT,
        reads = os.path.join(TOTAL_READS_DIR, "total_reads_{tw}.tsv.gz")
    output:
        tsv = os.path.join(OUTDIR, "negbinom_parameters_{tw}_all_and_by_population.tsv.gz"),
        rds = os.path.join(OUTDIR, "negbinom_models_{tw}_all_and_by_population.rds")
    log:
        os.path.join(OUTDIR, "logs", "{tw}.log")
    shell:
        """
        Rscript {input.script} "{wildcards.tw}" > {log} 2>&1
        """
