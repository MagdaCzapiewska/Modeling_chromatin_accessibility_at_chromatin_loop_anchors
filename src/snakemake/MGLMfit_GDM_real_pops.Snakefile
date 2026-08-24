import os

configfile: "config/config.yml"

RESULTSDIR = config["paths"]["resultsdir"]
DATADIR = config["paths"]["datadir"]

REAL_COUNTS_DIR = os.path.join(RESULTSDIR, "counts", "counts_in_anchors")
MGLM_REAL_POPS_BASE = os.path.join(RESULTSDIR, "MGLMfit_GDM", "real_data", "pops")
FIT_RSCRIPT = os.path.join(config['paths']['Rsrcdir'], "MGLMfit_GDM.R")

ALL_LOOPS = config["loops"].split(",")
CURRENT_TW = config["tw"]
CURRENT_POP = config["pop"]

INITS = ["default", "1e-6"]

rule all:
    input:
        [
            os.path.join(MGLM_REAL_POPS_BASE, f"init_{init}", f"fit_real_{CURRENT_TW}_{loop}_pop_{CURRENT_POP}.rds")
            for init in INITS
            for loop in ALL_LOOPS
        ]

rule fit_mglm_gdm_real_pop:
    input:
        script = FIT_RSCRIPT,
        counts = os.path.join(REAL_COUNTS_DIR, f"reads_{CURRENT_TW}_{{loop}}.tsv.gz")
    output:
        rds = os.path.join(MGLM_REAL_POPS_BASE, "init_{init}", f"fit_real_{CURRENT_TW}_{{loop}}_pop_{CURRENT_POP}.rds")
    params:
        custom_log = os.path.join(MGLM_REAL_POPS_BASE, "init_{init}", f"fit_GDM_{CURRENT_TW}_pop_{CURRENT_POP}.log")
    shell:
        """
        Rscript {input.script} \
            "{wildcards.init}" \
            "pop" \
            "{CURRENT_POP}" \
            "{input.counts}" \
            "{output.rds}" \
            "{params.custom_log}" > /dev/null 2>&1
        """
