import os

configfile: "config/config.yml"

RESULTSDIR = config["paths"]["resultsdir"]
# SNDIR = os.path.join(RESULTSDIR, "synthetic_data")
FIT_BASE_DIR = RESULTSDIR

OUT_BASE_DIR = os.path.join(RESULTSDIR, "MGLMfit_GDM_cor", "synthetic_data")
RSCRIPT = os.path.join(config['paths']['Rsrcdir'], "MGLMfit_GDM_cor", "cor_synthetic.R")

RHOS = [-0.6, 0.0, 0.6]
INITS = ["default", "1e-4", "1e-5", "1e-6", "1e-7", "1e-8"]

rule all:
    input:
        expand(
            os.path.join(OUT_BASE_DIR, "rho_{rho_str}", "init_{init}", "cor.tsv.gz"),
            rho_str=[f"{r:.1f}" for r in RHOS],
            init=INITS
        )

rule calculate_synthetic_cor:
    input:
        script = RSCRIPT,
        funcs = os.path.join(config['paths']['Rsrcdir'], "MGLMfit_GDM_cor", "correlation_functions.R"),
        fit_path = FIT_BASE_DIR
        # syn_path = SNDIR
    output:
        os.path.join(OUT_BASE_DIR, "rho_{rho_str}", "init_{init}", "cor.tsv.gz")
    log:
        os.path.join(OUT_BASE_DIR, "logs", "cor_rho_{rho_str}_init_{init}.log")
    shell:
        """
        Rscript {input.script} \
            "{wildcards.rho_str}" \
            "{wildcards.init}" \
            "{output}" \
            "{input.fit_path}" > {log} 2>&1
        """
    # shell:
    #     """
    #     Rscript {input.script} \
    #         "{wildcards.rho_str}" \
    #         "{wildcards.init}" \
    #         "{output}" \
    #         "{input.fit_path}" \
    #         "{input.syn_path}" > {log} 2>&1
    #     """
