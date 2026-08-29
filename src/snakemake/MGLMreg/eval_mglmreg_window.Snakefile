import os

configfile: "config/config.yml"

RESULTSDIR = config["paths"]["resultsdir"]

TIME_WINDOWS = [
    "00-02", "02-04", "04-06", "06-08", "08-10",
    "10-12", "12-14", "14-16", "16-18", "18-20"
]

OUT_DIR = os.path.join(RESULTSDIR, "MGLMreg", "eval_0h+_time_default")
RSCRIPT = os.path.join(config["paths"]["Rsrcdir"], "MGLMreg", "eval_mglmreg_window.R")

rule all:
    input:
        expand(
            os.path.join(OUT_DIR, "eval_cor_{tw}.tsv.gz"),
            tw=TIME_WINDOWS
        )

rule evaluate_mglmreg_window:
    input:
        script = RSCRIPT,
        fit_cor = os.path.join(RESULTSDIR, "MGLMfit_GDM_cor", "real_data", "all", "init_1e-6", "cor_{tw}.tsv.gz"),
        summary_reg = os.path.join(RESULTSDIR, "MGLMreg", "summary_GDM_0h+_time_default.tsv")
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
            "{output.tsv}" > {log} 2>&1
        """
