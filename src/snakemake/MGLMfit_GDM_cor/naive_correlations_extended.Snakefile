import os

configfile: "config/config.yml"

RESULTSDIR = config["paths"]["resultsdir"]

SNDIR = os.path.join(RESULTSDIR, "synthetic_data_extended")

OUT_BASE_DIR = os.path.join(RESULTSDIR, "naive_correlation", "synthetic_data_extended")
RSCRIPT = os.path.join(config['paths']['Rsrcdir'], "MGLMfit_GDM_cor", "naive_cor_synthetic_extended.R")

#RHOS = [-1.0, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
SEEDS = list(range(20))

RHOS = [-0.6, 0.0, 0.6]
#SEEDS = list(range(1000))

rule all:
    input:
        expand(
            os.path.join(OUT_BASE_DIR, "rho_{rho_str}", "naive_cor_seed{seed}.tsv.gz"),
            rho_str=[f"{r:.1f}" for r in RHOS],
            seed=SEEDS
        )

rule calculate_synthetic_cor:
    input:
        script = RSCRIPT,
        marker = os.path.join(SNDIR, "rho_{rho_str}", "sim_seed{seed}.done")
    output:
        os.path.join(OUT_BASE_DIR, "rho_{rho_str}", r"naive_cor_seed{seed,\d+}.tsv.gz")
    log:
        os.path.join(OUT_BASE_DIR, "logs", "cor_rho_{rho_str}_seed{seed}.log")
    params:
        rho_val = lambda wildcards: wildcards.rho_str
    shell:
        """
        Rscript {input.script} \
            "{params.rho_val}" \
            "{wildcards.seed}" \
            "{output}" \
            "{SNDIR}" > {log} 2>&1
        """
