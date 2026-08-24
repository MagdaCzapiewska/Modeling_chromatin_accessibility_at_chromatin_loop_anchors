import os
import pandas as pd

configfile: "config/config.yml"

RESULTSDIR = config["paths"]["resultsdir"]
DATADIR = config["paths"]["datadir"]
REAL_COUNTS_DIR = os.path.join(RESULTSDIR, "counts", "counts_in_anchors")
MGLM_REAL_BASE = os.path.join(RESULTSDIR, "MGLMfit_GDM", "real_data", "all")
FIT_RSCRIPT = os.path.join(config['paths']['Rsrcdir'], "MGLMfit_GDM.R")

################################################################################

TIME_WINDOWS = ["00-02", "02-04", "04-06", "06-08", "08-10", "10-12", "12-14", "14-16", "16-18", "18-20"]

INITS = ["default", "1e-6"]

loops_file = os.path.join(DATADIR, "long_and_short_range_loops_D_mel.tsv")
loops_table = pd.read_csv(loops_file, sep="\t")
loop_ids = loops_table["loop_id"].dropna().unique().tolist()

################################################################################

rule all:
    input:
        expand(
            os.path.join(MGLM_REAL_BASE, "init_{init}", "fit_real_{tw}_{loop}.rds"),
            init=INITS,
            tw=TIME_WINDOWS,
            loop=loop_ids
        )

rule fit_mglm_gdm_real_all:
    input:
        script = FIT_RSCRIPT,
        counts = os.path.join(REAL_COUNTS_DIR, "reads_{tw}_{loop}.tsv.gz")
    output:
        rds = os.path.join(MGLM_REAL_BASE, "init_{init}", "fit_real_{tw}_{loop}.rds")
    params:
        custom_log = os.path.join(MGLM_REAL_BASE, "init_{init}", "fit_GDM_{tw}.log")
    shell:
        """
        Rscript {input.script} \
            "{wildcards.init}" \
            "all" \
            "NA" \
            "{input.counts}" \
            "{output.rds}" \
            "{params.custom_log}" > /dev/null 2>&1
        """
