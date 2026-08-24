import os

configfile: "config/config.yml"

RESULTSDIR = config["paths"]["resultsdir"]
SNDIR = os.path.join(RESULTSDIR, "synthetic_data")
MGLM_BASE = os.path.join(RESULTSDIR, "MGLMfit_GDM", "synthetic_data")
FIT_RSCRIPT = os.path.join(config['paths']['Rsrcdir'], "MGLMfit_GDM.R")

RHOS = [-0.6, 0.0, 0.6]

N_CELLS = [
    500, 1000, 1500, 2000, 2500, 3000, 3500, 4000, 5000, 
    10000, 15000, 20000, 25000, 30000, 35000, 40000, 45000, 
    50000, 55000, 60000, 80000, 100000, 120000, 140000, 160000
]

ALPHA = [0.15]
BETA = [600]

MUS = [1000, 2000, 3000, 4000, 5000, 6000]

SIZE_NEGBINOMS = ["inf", "fixed", "0.5", "1.0", "1.5", "2.0", "2.5", "3.0", "3.5", "10"]

INITS = ["default", "1e-4", "1e-5", "1e-6", "1e-7", "1e-8"]

################################################################################

rule all:
    input:
        expand(
            os.path.join(MGLM_BASE, "rho_{rho_str}", "init_{init}", "fit_synthetic_n{n}_mu{mu}_sizeNB{size_nb}_alpha{alpha}_beta{beta}.rds"),
            rho_str=[f"{r:.1f}" for r in RHOS],
            init=INITS,
            n=N_CELLS,
            mu=MUS,
            size_nb=SIZE_NEGBINOMS,
            alpha=ALPHA,
            beta=BETA
        )

rule fit_mglm_gdm:
    input:
        script = FIT_RSCRIPT,
        counts = os.path.join(SNDIR, "rho_{rho_str}", r"synthetic_counts_n{n,\d+}_mu{mu,\d+}_sizeNB{size_nb}_alpha{alpha}_beta{beta}.tsv.gz")
    output:
        rds = os.path.join(MGLM_BASE, "rho_{rho_str}", "init_{init}", r"fit_synthetic_n{n,\d+}_mu{mu,\d+}_sizeNB{size_nb}_alpha{alpha}_beta{beta}.rds")
    params:
        custom_log = os.path.join(MGLM_BASE, "rho_{rho_str}", "init_{init}", "fit_GDM.log")
    shell:
        """
        Rscript {input.script} \
            "{wildcards.init}" \
            "all" \
            "NA" \
            "{input.counts}" \
            "{output.rds}" \
            "{params.custom_log}" > /dev/null 2>&1
        """
