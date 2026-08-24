import os

configfile: "config/config.yml"

RESULTSDIR = config["paths"]["resultsdir"]
INITS = ["default", "1e-6"]
TARGET_TW = ["06-08", "10-12", "14-16"]

DIR_REAL_ALL = os.path.join(RESULTSDIR, "MGLMfit_GDM_cor", "real_data", "all")
DIR_REAL_POPS = os.path.join(RESULTSDIR, "MGLMfit_GDM_cor", "real_data", "pops")

PLOT_SCRIPT = os.path.join(config['paths']['Rsrcdir'], "MGLMfit_GDM_cor", "plot_activity_correlations.R")

rule all:
    input:
        expand(os.path.join(DIR_REAL_ALL, "init_{init}", "plots_activity_report.pdf"), init=INITS),
        expand(os.path.join(DIR_REAL_POPS, "init_{init}", "plots_activity_report.pdf"), init=INITS)

rule plot_activity_all_correlations:
    input:
        script = PLOT_SCRIPT,
        files = expand(os.path.join(DIR_REAL_ALL, "init_{{init}}", "cor_{tw}.tsv.gz"), tw=TARGET_TW)
    output:
        os.path.join(DIR_REAL_ALL, "init_{init}", "plots_activity_report.pdf")
    log:
        os.path.join(DIR_REAL_ALL, "logs", "plotting_activity_init_{init}.log")
    shell:
        """
        INPUT_DIR=$(dirname "{input.files[0]}")
        
        Rscript {input.script} \
            "all" \
            "{wildcards.init}" \
            "$INPUT_DIR" \
            "{output}" > {log} 2>&1
        """

rule plot_activity_pops_correlations:
    input:
        script = PLOT_SCRIPT,
        files = expand(os.path.join(DIR_REAL_POPS, "init_{{init}}", "cor_{tw}.tsv.gz"), tw=TARGET_TW)
    output:
        os.path.join(DIR_REAL_POPS, "init_{init}", "plots_activity_report.pdf")
    log:
        os.path.join(DIR_REAL_POPS, "logs", "plotting_activity_init_{init}.log")
    shell:
        """
        INPUT_DIR=$(dirname "{input.files[0]}")
        
        Rscript {input.script} \
            "pops" \
            "{wildcards.init}" \
            "$INPUT_DIR" \
            "{output}" > {log} 2>&1
        """
