import os

configfile: "config/config.yml"

RESULTSDIR = config["paths"]["resultsdir"]
IN_COR_DIR = os.path.join(RESULTSDIR, "naive_correlation", "synthetic_data_extended")
OUT_PLOT_DIR = os.path.join(RESULTSDIR, "plots_extended", "naive_correlations")
PLOT_SCRIPT = os.path.join(config['paths']['Rsrcdir'], "MGLMfit_GDM_cor", "plot_synthetic_extended_naive.R")

RHOS = [-1.0, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
SEEDS = list(range(20))

rule all:
    input:
        expand(
            os.path.join(OUT_PLOT_DIR, "rho_{rho_str}", "{pdf_type}.pdf"),
            rho_str=[f"{r:.1f}" for r in RHOS],
            pdf_type=["variable_N", "variable_MU", "variable_SIZE"]
        )

rule plot_extended_naive_distributions:
    input:
        script = PLOT_SCRIPT,
        cor_files = expand(os.path.join(IN_COR_DIR, "rho_{rho_str}", "naive_cor_seed{seed}.tsv.gz"), seed=SEEDS, allow_missing=True)
    output:
        pdf_n = os.path.join(OUT_PLOT_DIR, "rho_{rho_str}", "variable_N.pdf"),
        pdf_mu = os.path.join(OUT_PLOT_DIR, "rho_{rho_str}", "variable_MU.pdf"),
        pdf_sz = os.path.join(OUT_PLOT_DIR, "rho_{rho_str}", "variable_SIZE.pdf")
    log:
        os.path.join(OUT_PLOT_DIR, "logs", "plotting_rho_{rho_str}.log")
    params:
        in_dir = IN_COR_DIR,
        out_dir = lambda wildcards: os.path.join(OUT_PLOT_DIR, f"rho_{wildcards.rho_str}"),
        rho_val = lambda wildcards: wildcards.rho_str
    shell:
        """
        Rscript {input.script} \
            "{params.in_dir}" \
            "{params.rho_val}" \
            "{params.out_dir}" > {log} 2>&1
        """
