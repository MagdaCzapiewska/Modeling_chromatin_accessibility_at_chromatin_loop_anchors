import os

configfile: "config/config.yml"

RESULTSDIR = config["paths"]["resultsdir"]

SNDIR = os.path.join(RESULTSDIR, "synthetic_data")

FITDIR = os.path.join(RESULTSDIR, "MGLMfit_init_1e-6")

OUTDIR = os.path.join(RESULTSDIR, "MGLMfit_init_1e-6_check_and_cor")

RSCRIPT = os.path.join(config['paths']['Rsrcdir'], "MGLMfit_GDM_check_and_corr.R")

RHOS = [-0.6, 0.0, 0.6]

rule all:
    input:
        expand(
            os.path.join(OUTDIR, "check_and_cor_rho_{rho_str}.tsv.gz"),
            rho_str=[f"{r:.1f}" for r in RHOS]
        )

rule calculate_correlation_grid:
    input:
        fit_path = FITDIR,
        syn_path = SNDIR
    output:
        os.path.join(OUTDIR, "check_and_cor_rho_{rho_str}.tsv.gz")
    log:
        os.path.join(OUTDIR, "logs", "cor_rho_{rho_str}.log")
    params:
        rho_val = lambda wildcards: wildcards.rho_str
    shell:
        """
        Rscript {RSCRIPT} \
            {params.rho_val} \
            {output} \
            {input.fit_path} \
            {input.syn_path} > {log} 2>&1
        """
