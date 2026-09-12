import os

configfile: "config/config.yml"

RESULTSDIR = config["paths"]["resultsdir"]
SNDIR = os.path.join(RESULTSDIR, "synthetic_data_extended_plus")
MGLM_BASE = os.path.join(RESULTSDIR, "MGLMfit_GDM", "synthetic_data_extended_plus")
FIT_RSCRIPT = os.path.join(config['paths']['Rsrcdir'], "MGLMfit_GDM_extended_plus.R")

RHOS = [-0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6]
SEEDS = list(range(1000))
INITS = ["1e-6"]

rule all:
    input:
        expand(
            os.path.join(MGLM_BASE, "rho_{rho_str}", "init_{init}", "fit_seed{seed}.done"),
            rho_str=[f"{r:.1f}" for r in RHOS],
            init=INITS,
            seed=SEEDS
        )

rule fit_mglm_gdm_extended:
    input:
        synthetic_file = os.path.join(SNDIR, "rho_{rho_str}", r"synthetic_counts_seed{seed,\d+}.tsv.gz")
    output:
        done = os.path.join(MGLM_BASE, "rho_{rho_str}", "init_{init}", r"fit_seed{seed,\d+}.done")
    log:
        os.path.join(MGLM_BASE, "logs", "rho_{rho_str}", r"fit_rho_{rho_str}_init_{init}_seed{seed}.log")
    params:
        custom_log = os.path.join(MGLM_BASE, "rho_{rho_str}", "init_{init}", "fit_GDM_shared.log"),
        rho_val = lambda wildcards: wildcards.rho_str,
        script = FIT_RSCRIPT
    shell:
        """
        export OMP_NUM_THREADS=1
        export MKL_NUM_THREADS=1
        export OPENBLAS_NUM_THREADS=1
        export VECLIB_MAXIMUM_THREADS=1
        export NUMEXPR_NUM_THREADS=1

        export MPI_NUM_THREADS=1
        export OMPI_NUM_THREADS=1
        
        Rscript {params.script} \
            "{wildcards.init}" \
            "{params.rho_val}" \
            "{wildcards.seed}" \
            "{SNDIR}" \
            "{MGLM_BASE}" \
            "{params.custom_log}" > {log} 2>&1
        """
