import os

configfile: "config/config.yml"

RESULTSDIR = config["paths"]["resultsdir"]

DIR_NAIVE_REAL_ALL  = os.path.join(RESULTSDIR, "naive_correlations", "real_data", "all")
DIR_NAIVE_REAL_POPS = os.path.join(RESULTSDIR, "naive_correlations", "real_data", "pops")
DIR_NAIVE_SYN       = os.path.join(RESULTSDIR, "naive_correlations", "synthetic_data")

SCRIPT_REAL = os.path.join(config['paths']['Rsrcdir'], "MGLMfit_GDM_cor", "plot_naive_real.R")
SCRIPT_SYN  = os.path.join(config['paths']['Rsrcdir'], "MGLMfit_GDM_cor", "plot_naive_synthetic.R")

TIME_WINDOWS = ["00-02", "02-04", "04-06", "06-08", "08-10", "10-12", "12-14", "14-16", "16-18", "18-20"]
RHOS = [-0.6, 0.0, 0.6]

rule all:
    input:
        os.path.join(DIR_NAIVE_REAL_ALL, "plots_report_naive.pdf"),
        os.path.join(DIR_NAIVE_REAL_POPS, "plots_report_naive.pdf"),
        expand(os.path.join(DIR_NAIVE_SYN, "rho_{rho_str}", "plots_report_naive.pdf"), rho_str=[f"{r:.1f}" for r in RHOS])

rule plot_naive_real_all_correlations:
    input:
        script = SCRIPT_REAL,
        files = expand(os.path.join(DIR_NAIVE_REAL_ALL, "cor_{tw}.tsv.gz"), tw=TIME_WINDOWS)
    output:
        os.path.join(DIR_NAIVE_REAL_ALL, "plots_report_naive.pdf")
    log:
        os.path.join(DIR_NAIVE_REAL_ALL, "logs", "plotting_naive.log")
    shell:
        """
        INPUT_DIR=$(dirname "{input.files[0]}")
        Rscript {input.script} "all" "$INPUT_DIR" "{output}" > {log} 2>&1
        """

rule plot_naive_real_pops_correlations:
    input:
        script = SCRIPT_REAL,
        files = expand(os.path.join(DIR_NAIVE_REAL_POPS, "cor_{tw}.tsv.gz"), tw=TIME_WINDOWS)
    output:
        os.path.join(DIR_NAIVE_REAL_POPS, "plots_report_naive.pdf")
    log:
        os.path.join(DIR_NAIVE_REAL_POPS, "logs", "plotting_naive.log")
    shell:
        """
        INPUT_DIR=$(dirname "{input.files[0]}")
        Rscript {input.script} "pops" "$INPUT_DIR" "{output}" > {log} 2>&1
        """

rule plot_naive_synthetic_distribution:
    input:
        script = SCRIPT_SYN,
        tsv = os.path.join(DIR_NAIVE_SYN, "rho_{rho_str}", "cor.tsv.gz")
    output:
        os.path.join(DIR_NAIVE_SYN, "rho_{rho_str}", "plots_report_naive.pdf")
    log:
        os.path.join(DIR_NAIVE_SYN, "logs", "plotting_rho_{rho_str}.log")
    shell:
        """
        Rscript {input.script} "{input.tsv}" "{wildcards.rho_str}" "{output}" > {log} 2>&1
        """
