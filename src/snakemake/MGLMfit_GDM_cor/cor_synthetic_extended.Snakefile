import os

configfile: "config/config.yml"

RESULTSDIR = config["paths"]["resultsdir"]
SNDIR = os.path.join(RESULTSDIR, "synthetic_data_extended")
FIT_BASE_DIR = RESULTSDIR

OUT_BASE_DIR = os.path.join(RESULTSDIR, "MGLMfit_GDM_cor", "synthetic_data_extended")
RSCRIPT = os.path.join(config['paths']['Rsrcdir'], "MGLMfit_GDM_cor", "cor_synthetic_extended.R")

#RHOS = [-1.0, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
#SEEDS = list(range(20))

RHOS = [-0.6, 0.0, 0.6]
SEEDS = list(range(20, 1000))

INITS = ["1e-6"]

rule all:
    input:
        expand(
            os.path.join(OUT_BASE_DIR, "rho_{rho_str}", "init_{init}", "cor_seed{seed}.tsv.gz"),
            rho_str=[f"{r:.1f}" for r in RHOS],
            init=INITS,
            seed=SEEDS
        )

rule calculate_synthetic_cor_extended:
    input:
        #script = RSCRIPT,
        funcs = os.path.join(config['paths']['Rsrcdir'], "MGLMfit_GDM_cor", "correlation_functions.R"),
        fit_marker = os.path.join(RESULTSDIR, "MGLMfit_GDM", "synthetic_data_extended", "rho_{rho_str}", "init_{init}", "fit_seed{seed}.done")
    output:
        tsv = os.path.join(OUT_BASE_DIR, "rho_{rho_str}", "init_{init}", r"cor_seed{seed,\d+}.tsv.gz")
    log:
        os.path.join(OUT_BASE_DIR, "logs", "rho_{rho_str}", "cor_rho_{rho_str}_init_{init}_seed{seed}.log")
    params:
        fit_path = FIT_BASE_DIR,
        syn_path = SNDIR,
        rho_val = lambda wildcards: wildcards.rho_str,
        script = RSCRIPT
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
            "{params.rho_val}" \
            "{wildcards.init}" \
            "{wildcards.seed}" \
            "{output.tsv}" \
            "{params.fit_path}" \
            "{params.syn_path}" > {log} 2>&1
        """
