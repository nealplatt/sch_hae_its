#!/bin/bash

# Activate conda environment
conda activate mito_assembly

# Define directories and commands
PROJ_DIR="/master/nplatt/sch_hae_its-nigeria"
RESULTS_DIR="$PROJ_DIR/results"
QSUB="qsub -V -cwd -S /bin/bash -q all.q -j y -p -1023"
CONDA="conda activate $PROJ_DIR/envs/mito_assembly"
SEED="${RESULTS_DIR}/mito_read_count/reference_mitochondria.fasta"

# Change to the appropriate directory
cd "$RESULTS_DIR/mito_assembly"

# Add animal_mt to the database (if not already added)
get_organelle_config.py --add animal_mt

# Ensure logs directory exists
mkdir -p "$RESULTS_DIR/mito_assembly/logs" "$RESULTS_DIR/mito_assembly/get_organelle"

# Loop through samples in samples.list
while read -r SAMPLE; do
    
    # Define input file paths
    FQ="$RESULTS_DIR/mito_read_count/${SAMPLE}.ambig-all.mapped.fq.gz"

    OUT_DIR="$RESULTS_DIR/mito_assembly/get_organelle/$SAMPLE"
    LOG_FILE="$RESULTS_DIR/mito_assembly/logs/$SAMPLE.log"
    
    # Construct the command
    CMD="conda run -n mito_assembly get_organelle_from_reads.py \
        -u $FQ \
        -t 12 \
        -R 10 \
        -s $SEED \
        -k 21,45,65,85 \
        -F animal_mt \
        -o $OUT_DIR"

    # Submit the job
    echo "$CMD" | $QSUB -pe smp 12 -N mito_$SAMPLE -o "$LOG_FILE"
done < samples.list

