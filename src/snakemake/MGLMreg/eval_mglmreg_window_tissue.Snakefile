import os

configfile: "config/config.yml"

RESULTSDIR = config["paths"]["resultsdir"]

TIME_WINDOWS = [
    "10-12", "12-14", "14-16", "16-18", "18-20"
]

OUT_DIR = os.path.join(RESULTSDIR, "MGLMreg", "eval_10h+_time_tissue_default")
RSCRIPT = os.path.join(config["paths"]["Rsrcdir"], "MGLMreg", "eval_mglmreg_window_tissue.R")

rule all:
    input:
        expand(
            os.path.join(OUT_DIR, "eval_cor_{tw}.tsv.gz"),
            tw=TIME_WINDOWS
        )

rule evaluate_mglmreg_window_tissue:
    input:
        script = RSCRIPT,
        pop_map = os.path.join(RESULTSDIR, "cluster_pop_mapping.tsv.gz"),
        fit_cor = os.path.join(RESULTSDIR, "MGLMfit_GDM_cor", "real_data", "pops", "init_1e-6", "cor_{tw}.tsv.gz")
    output:
        tsv = os.path.join(OUT_DIR, "eval_cor_{tw}.tsv.gz")
    log:
        os.path.join(OUT_DIR, "logs", "eval_{tw}.log")
    threads: 1
    shell:
        """
        export OMP_NUM_THREADS=1
        export OPENBLAS_NUM_THREADS=1

        Rscript {input.script} \
            "{wildcards.tw}" \
            "{output.tsv}" \
            "default" > {log} 2>&1
        """
