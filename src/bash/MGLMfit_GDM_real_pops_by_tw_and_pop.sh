#!/bin/bash

LOOPS_FILE="data/long_and_short_range_loops_D_mel.tsv"
CARDINALITY_FILE="results/cluster_cardinality.tsv.gz"
SNAKEFILE="src/snakemake/MGLMfit_GDM_real_pops.Snakefile"
CORES=32

echo "=== Parsing loop IDs ==="
mapfile -t ALL_LOOPS < <(awk -F'\t' 'NR==1 {for(i=1;i<=NF;i++) if($i=="loop_id") col=i; next} col && $col!="" && !seen[$col]++ {print $col}' "$LOOPS_FILE")
TOTAL_LOOPS=${#ALL_LOOPS[@]}

ALL_LOOPS_STRING=$(IFS=,; echo "${ALL_LOOPS[*]}")

echo "=== Parsing population numbers ==="
mapfile -t VALID_COMBINATIONS < <(zcat "$CARDINALITY_FILE" | awk -F'\t' '
    NR==1 {
        for(i=2;i<=NF;i++) tw[i]=$i
        next
    }
    {
        pop=$1
        for(i=2;i<=NF;i++) {
            if($i > 0) print tw[i]","pop
        }
    }
')

TOTAL_COMBOS=${#VALID_COMBINATIONS[@]}
echo "Number of loops found: $TOTAL_LOOPS"
echo "Number of combinations (time window, population) found: $TOTAL_COMBOS"
echo "--------------------------------------------------"

COMBO_COUNTER=0

for combo in "${VALID_COMBINATIONS[@]}"; do
    IFS=',' read -r current_tw current_pop <<< "$combo"
    ((COMBO_COUNTER++))
    
    echo "=========================================================================="
    echo " Step [$COMBO_COUNTER/$TOTAL_COMBOS] -> TW: $current_tw | POPULATION: $current_pop"
    echo "=========================================================================="
    
    echo "[$(date +'%Y-%m-%d %H:%M:%S')] Starting Snakemake for all $TOTAL_LOOPS loops..."

    snakemake \
        -s "$SNAKEFILE" \
        --cores "$CORES" \
        --config loops="$ALL_LOOPS_STRING" tw="$current_tw" pop="$current_pop"
        
    echo "[$(date +'%Y-%m-%d %H:%M:%S')] Finished step for TW $current_tw | POP $current_pop"
done

echo "=== Pipeline finished without error ==="
