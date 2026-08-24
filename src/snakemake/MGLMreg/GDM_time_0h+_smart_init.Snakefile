import os
import pandas as pd

configfile: "config/config.yml"

RESULTSDIR = config["paths"]["resultsdir"]
DATADIR = os.path.join(RESULTSDIR, "counts", "counts_in_anchors")
OUTPUTDIR = os.path.join(RESULTSDIR, "MGLMreg", "GDM_time_0h+")
RSCRIPT = os.path.join(config["paths"]["Rsrcdir"], "MGLMreg", "GDM_time_0h+_smart_init.R")
LOOPS_FILE = os.path.join(config["paths"]["datadir"], "long_and_short_range_loops_D_mel.tsv")

loops_df = pd.read_csv(LOOPS_FILE, sep="\t")
all_loops = sorted(loops_df["loop_id"].unique().tolist())

n_take = 1000
use_all = True
time_windows = ["00-02","02-04","04-06","06-08","08-10",
                "10-12","12-14","14-16","16-18","18-20"]

def out_file(loop, n_take, use_all):
    if use_all:
        return os.path.join(OUTPUTDIR, f"fit_reg_{loop}.rds")
    else:
        return os.path.join(OUTPUTDIR, f"fit_reg_{loop}_n_points_{int(n_take)}.rds")

# ---- rule all ----
rule all:
    input:
        expand(out_file("{loop}", n_take, use_all), loop=all_loops)

# ---- rule: fit_gdm ----
rule fit_gdm:
    input:
        script = RSCRIPT,
        data = expand(os.path.join(DATADIR, "reads_{tw}_{{loop}}.tsv.gz"), tw=time_windows)
    output:
        rds = out_file("{loop}", n_take, use_all)
    params:
        loop = "{loop}",
        all_flag = str(use_all).upper(),
        n_take = n_take
    threads: 1
    shell:
        """
        OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
        Rscript {input.script} {params.loop} {output.rds} {params.all_flag} {params.n_take}
        """
