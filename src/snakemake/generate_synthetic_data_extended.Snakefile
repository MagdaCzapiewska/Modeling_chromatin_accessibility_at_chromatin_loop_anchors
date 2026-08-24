import os

configfile: "config/config.yml"

RESULTSDIR = config["paths"]["resultsdir"]
OUTDIR = os.path.join(RESULTSDIR, "synthetic_data_extended")
RSCRIPT = os.path.join(config['paths']['Rsrcdir'], "generate_synthetic_data_extended.R")

#RHOS = [-1.0, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
#SEEDS = list(range(20))

RHOS = [-0.6, 0.0, 0.6]
SEEDS = list(range(1000))

rule all:
    input:
        expand(
            os.path.join(OUTDIR, "rho_{rho_str}", "sim_seed{seed}.done"),
            rho_str=[f"{r:.1f}" for r in RHOS],
            seed=SEEDS
        )

rule generate_synthetic_data:
    output:
        os.path.join(OUTDIR, "rho_{rho_str}", r"sim_seed{seed,\d+}.done")
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
