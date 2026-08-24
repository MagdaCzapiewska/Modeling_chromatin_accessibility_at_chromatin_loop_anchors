import os
import glob

# Load the project configuration
configfile: "config/config.yml"

RESULTSDIR = config["paths"]["resultsdir"]
SCRIPT_TRAJECTORY = os.path.join(config['paths']['Rsrcdir'], "MGLMfit_GDM_cor", "plot_correlation_trajectories.R")

TIME_WINDOWS = ["00-02", "02-04", "04-06", "06-08", "08-10", "10-12", "12-14", "14-16", "16-18", "18-20"]

# ------------------------------------------------------------------------------
# DYNAMIC INITS DISCOVERY
# ------------------------------------------------------------------------------
# Dynamically find all initialization subdirectories present in the GDM results path
GDM_BASE_PATH = os.path.join(RESULTSDIR, "MGLMfit_GDM_cor", "real_data", "all")

if os.path.exists(GDM_BASE_PATH):
    # Captures 'default', 'random', etc., from 'real_data/all/init_*/'
    INIT_DIRS = [
        os.path.basename(d).replace("init_", "") 
        for d in glob.glob(os.path.join(GDM_BASE_PATH, "init_*")) 
        if os.path.isdir(d)
    ]
else:
    INIT_DIRS = ["default"] # Fallback if directories aren't built yet

# If no directories were matched but the path exists, ensure at least 'default' is targeted
if not INIT_DIRS:
    INIT_DIRS = ["default"]

# ------------------------------------------------------------------------------
# TARGETS RULE
# ------------------------------------------------------------------------------
rule all_trajectories:
    input:
        # Target for the Naive Baseline data
        os.path.join(RESULTSDIR, "naive_correlations", "real_data", "all", "trajectories_report.pdf"),
        # Dynamic targets for all discovered GDM initialization methods
        expand(os.path.join(GDM_BASE_PATH, "init_{init}", "trajectories_report.pdf"), init=INIT_DIRS)


# ------------------------------------------------------------------------------
# WORKFLOW RULES
# ------------------------------------------------------------------------------

rule plot_naive_all_trajectories:
    """
    Generates a trajectory density cloud report for ALL loops using naive raw counts.
    """
    input:
        script = SCRIPT_TRAJECTORY,
        files = expand(os.path.join(RESULTSDIR, "naive_correlations", "real_data", "all", "cor_{tw}.tsv.gz"), tw=TIME_WINDOWS)
    output:
        pdf = os.path.join(RESULTSDIR, "naive_correlations", "real_data", "all", "trajectories_report.pdf")
    log:
        log = os.path.join(RESULTSDIR, "naive_correlations", "real_data", "all", "logs", "trajectories_generation.log")
    shell:
        """
        INPUT_DIR=$(dirname "{input.files[0]}")
        Rscript {input.script} "naive" "$INPUT_DIR" "{output.pdf}" > {log.log} 2>&1
        """


rule plot_gdm_all_trajectories:
    """
    Generates a trajectory density cloud report for ALL loops across any GDM initialization model.
    """
    input:
        script = SCRIPT_TRAJECTORY,
        files = expand(os.path.join(GDM_BASE_PATH, "init_{init}", "cor_{tw}.tsv.gz"), tw=TIME_WINDOWS, allow_missing=True)
    output:
        pdf = os.path.join(GDM_BASE_PATH, "init_{init}", "trajectories_report.pdf")
    log:
        log = os.path.join(GDM_BASE_PATH, "init_{init}", "logs", "trajectories_generation.log")
    shell:
        """
        INPUT_DIR=$(dirname "{input.files[0]}")
        Rscript {input.script} "gdm" "$INPUT_DIR" "{output.pdf}" > {log.log} 2>&1
        """
