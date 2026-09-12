import os

configfile: "config/config.yml"

RESULTSDIR = config["paths"]["resultsdir"]

OUT_BASE_DIR = os.path.join(RESULTSDIR, "MGLMfit_GDM_cor", "synthetic_data")
PLOT_SCRIPT = os.path.join(config['paths']['Rsrcdir'], "MGLMfit_GDM_cor", "plot_synthetic_correlations.R")

RHOS = [-0.6, 0.0, 0.6]
INITS = ["default", "1e-4", "1e-5", "1e-6", "1e-7", "1e-8"]

rule all:
    input:
        expand(
            os.path.join(OUT_BASE_DIR, "rho_{rho_str}", "init_{init}", "plots_report.pdf"),
            rho_str=[f"{r:.1f}" for r in RHOS],
            init=INITS
        )

rule plot_synthetic_distribution:
    input:
        script = PLOT_SCRIPT,
        tsv = os.path.join(OUT_BASE_DIR, "rho_{rho_str}", "init_{init}", "cor.tsv.gz")
    output:
        os.path.join(OUT_BASE_DIR, "rho_{rho_str}", "init_{init}", "plots_report.pdf")
    log:
        os.path.join(OUT_BASE_DIR, "logs", "plotting_rho_{rho_str}_init_{init}.log")
    shell:
        """
        Rscript {input.script} \
            "{input.tsv}" \
            "{wildcards.rho_str}" \
            "{output}" > {log} 2>&1
        """
