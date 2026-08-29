import os

configfile: "config/config.yml"

RESULTSDIR = config["paths"]["resultsdir"]
OUTDIR = os.path.join(RESULTSDIR, "synthetic_data_extended_plus")
RSCRIPT = os.path.join(config['paths']['Rsrcdir'], "generate_synthetic_data_extended_plus.R")

RHOS = [-0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6]
SEEDS = list(range(1000))

rule all:
    input:
        expand(
            os.path.join(OUTDIR, "rho_{rho_str}", "synthetic_counts_seed{seed}.tsv.gz"),
            rho_str=[f"{r:.1f}" for r in RHOS],
            seed=SEEDS
        )

rule generate_synthetic_data:
    output:
        os.path.join(OUTDIR, "rho_{rho_str}", r"synthetic_counts_seed{seed,\d+}.tsv.gz")
    log:
        os.path.join(OUTDIR, "logs", "rho_{rho_str}", "sim_seed{seed}.log")
    params:
        rho_val = lambda wildcards: wildcards.rho_str
    shell:
        """
        Rscript {RSCRIPT} \
            {params.rho_val} \
            {wildcards.seed} > {log} 2>&1
        """
