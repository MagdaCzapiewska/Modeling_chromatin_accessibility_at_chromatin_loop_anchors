import os

configfile: "config/config.yml"

RESULTSDIR = config["paths"]["resultsdir"]
# REAL_COUNTS_DIR = os.path.join(RESULTSDIR, "counts", "counts_in_anchors")
FIT_BASE_DIR = RESULTSDIR

OUT_REAL_POPS = os.path.join(RESULTSDIR, "MGLMfit_GDM_cor", "real_data", "pops")
RSCRIPT = os.path.join(config['paths']['Rsrcdir'], "MGLMfit_GDM_cor", "cor_real_pops.R")

TIME_WINDOWS = ["00-02", "02-04", "04-06", "06-08", "08-10", "10-12", "12-14", "14-16", "16-18", "18-20"]
INITS = ["default", "1e-6"]

rule all:
    input:
        expand(
            os.path.join(OUT_REAL_POPS, "init_{init}", "cor_{tw}.tsv.gz"),
            init=INITS,
            tw=TIME_WINDOWS
        )

rule calculate_real_pops_cor:
    input:
        script = RSCRIPT,
        funcs = os.path.join(config['paths']['Rsrcdir'], "MGLMfit_GDM_cor", "correlation_functions.R"),
        fit_path = FIT_BASE_DIR,
        # counts = REAL_COUNTS_DIR
        pop_map = os.path.join(RESULTSDIR, "cluster_pop_mapping.tsv.gz"),
        cardinality = os.path.join(RESULTSDIR, "cluster_cardinality.tsv.gz")
    output:
        os.path.join(OUT_REAL_POPS, "init_{init}", "cor_{tw}.tsv.gz")
    log:
        os.path.join(OUT_REAL_POPS, "logs", "cor_{tw}_init_{init}.log")
    shell:
        """
        Rscript {input.script} \
            "{wildcards.tw}" \
            "{wildcards.init}" \
            "{output}" \
            "{input.fit_path}" > {log} 2>&1
        """
    # shell:
    #     """
    #     Rscript {input.script} \
    #         "{wildcards.tw}" \
    #         "{wildcards.init}" \
    #         "{output}" \
    #         "{input.fit_path}" \
    #         "{input.counts}" > {log} 2>&1
    #     """
