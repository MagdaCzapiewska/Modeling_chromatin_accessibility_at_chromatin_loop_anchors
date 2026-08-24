import os

configfile: "config/config.yml"

RESULTSDIR = config["paths"]["resultsdir"]

OUTDIR = os.path.join(RESULTSDIR, "synthetic_data")
RSCRIPT = os.path.join(config['paths']['Rsrcdir'], "generate_synthetic_data.R")

RHOS = [-0.6, 0.0, 0.6]

N_CELLS = [
    500, 1000, 1500, 2000, 2500, 3000, 3500, 4000, 5000, 
    10000, 15000, 20000, 25000, 30000, 35000, 40000, 45000, 
    50000, 55000, 60000, 80000, 100000, 120000, 140000, 160000
]

ALPHA = [0.15]
BETA = [600]

MUS = [1000, 2000, 3000, 4000, 5000, 6000]

# "fixed" = all total reads equal to mu, "inf" = Poisson
SIZE_NEGBINOMS = ["inf", "fixed", "0.5", "1.0", "1.5", "2.0", "2.5", "3.0", "3.5", "10"]

def get_rho_str(wildcards):
    """Formatting rho value to X.X (-0.6, 0.0)"""
    val = float(wildcards.rho)
    return f"{val:.1f}"

rule all:
    input:
        expand(
            os.path.join(OUTDIR, "rho_{rho_str}", "synthetic_counts_n{n}_mu{mu}_sizeNB{size_nb}_alpha{alpha}_beta{beta}.tsv.gz"),
            rho_str=[f"{r:.1f}" for r in RHOS],
            n=N_CELLS,
            mu=MUS,
            size_nb=SIZE_NEGBINOMS,
            alpha=ALPHA,
            beta=BETA
        )

rule generate_synthetic_data:
    output:
        os.path.join(OUTDIR, "rho_{rho_str}", "synthetic_counts_n{n,\d+}_mu{mu,\d+}_sizeNB{size_nb}_alpha{alpha}_beta{beta}.tsv.gz")
    log:
        os.path.join(OUTDIR, "logs", "rho_{rho_str}", "sim_n{n}_mu{mu}_sizeNB{size_nb}_alpha{alpha}_beta{beta}.log")
    params:
        # Mapping 'rho_str' to value used as argument for R script (for example -0.6)
        rho_val = lambda wildcards: wildcards.rho_str
    shell:
        """
        Rscript {RSCRIPT} \
            {params.rho_val} \
            {wildcards.n} \
            {wildcards.alpha} \
            {wildcards.beta} \
            {wildcards.mu} \
            {wildcards.size_nb} > {log} 2>&1
        """
