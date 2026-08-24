import os

configfile: "config/config.yml"

RESULTSDIR = config["paths"]["resultsdir"]
REAL_COUNTS_DIR = os.path.join(RESULTSDIR, "counts", "counts_in_anchors")
SNDIR = os.path.join(RESULTSDIR, "synthetic_data")

OUT_BASE = os.path.join(RESULTSDIR, "naive_correlations")

R_SCRIPT_REAL = os.path.join(config['paths']['Rsrcdir'], "MGLMfit_GDM_cor", "naive_cor_real.R")
R_SCRIPT_SYN  = os.path.join(config['paths']['Rsrcdir'], "MGLMfit_GDM_cor", "naive_cor_synthetic.R")

TIME_WINDOWS = ["00-02", "02-04", "04-06", "06-08", "08-10", "10-12", "12-14", "14-16", "16-18", "18-20"]
RHOS = [-0.6, 0.0, 0.6]

rule all:
    input:
        # 1. Real Data - All Cells
        expand(os.path.join(OUT_BASE, "real_data", "all", "cor_{tw}.tsv.gz"), tw=TIME_WINDOWS),
        # 2. Real Data - Populations
        expand(os.path.join(OUT_BASE, "real_data", "pops", "cor_{tw}.tsv.gz"), tw=TIME_WINDOWS),
        # 3. Synthetic Data
        expand(os.path.join(OUT_BASE, "synthetic_data", "rho_{rho_str}", "cor.tsv.gz"), rho_str=[f"{r:.1f}" for r in RHOS])

rule calculate_naive_real_all:
    input:
        script = R_SCRIPT_REAL,
        counts = REAL_COUNTS_DIR
    output:
        os.path.join(OUT_BASE, "real_data", "all", "cor_{tw}.tsv.gz")
    log:
        os.path.join(OUT_BASE, "logs", "naive_real_all_{tw}.log")
    shell:
        """
        Rscript {input.script} \
            "{wildcards.tw}" \
            "all" \
            "{output}" \
            "{input.counts}" > {log} 2>&1
        """

rule calculate_naive_real_pops:
    input:
        script = R_SCRIPT_REAL,
        counts = REAL_COUNTS_DIR,
        results_path = RESULTSDIR # Needed for cluster_cardinality.tsv.gz
    output:
        os.path.join(OUT_BASE, "real_data", "pops", "cor_{tw}.tsv.gz")
    log:
        os.path.join(OUT_BASE, "logs", "naive_real_pops_{tw}.log")
    shell:
        """
        Rscript {input.script} \
            "{wildcards.tw}" \
            "pops" \
            "{output}" \
            "{input.counts}" \
            "{input.results_path}" > {log} 2>&1
        """

rule calculate_naive_synthetic:
    input:
        script = R_SCRIPT_SYN,
        syn_path = SNDIR
    output:
        os.path.join(OUT_BASE, "synthetic_data", "rho_{rho_str}", "cor.tsv.gz")
    log:
        os.path.join(OUT_BASE, "logs", "naive_synthetic_rho_{rho_str}.log")
    shell:
        """
        Rscript {input.script} \
            "{wildcards.rho_str}" \
            "{output}" \
            "{input.syn_path}" > {log} 2>&1
        """
