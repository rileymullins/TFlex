#!/bin/bash
# step_3_run_span_peak_caller.sh
# ==============================================================================
# Configuration - Adjust these variables
# ==============================================================================

# --- Paths ---
# Path to the SPAN jar file
SPAN_JAR="/Users/rileymullins/Documents/span-2.0.6652.jar"

# Directory containing your input .qbed files
INPUT_DIR="/Users/rileymullins/Documents/test_scripts_for_multiplex_cc_analysis/TESTING_06302025/testing_07032025/full_dataset/full_dataset_insertion_maps/raw_unique_insertion_count_per_group_qbed"

# Directory where the output peak files will be saved
OUTPUT_DIR="/Users/rileymullins/Documents/span_analysis_results"

# Path to your chromosome sizes file
CHROM_SIZES="/Users/rileymullins/Documents/test_scripts_for_multiplex_cc_analysis/TESTING_06302025/testing_07032025/hg38.chrom.sizes"


# --- SPAN Parameters ---
# Memory to allocate to Java
JAVA_MEM="-Xmx16G"

# Bin size for the analysis
BIN_SIZE=50

# False Discovery Rate (FDR) cutoff
FDR="0.01"

# --- Concurrency Settings ---
# Set the maximum number of parallel jobs
MAX_JOBS=2


# ==============================================================================
# Script Logic - No need to edit below this line
# ==============================================================================

# Create the output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Initialize a counter for the jobs
job_count=0

# Find all .qbed files in the input directory and loop through them
for QBED_FILE in "$INPUT_DIR"/*.qbed; do
    # Get the base filename without the path and extension
    BASENAME=$(basename "$QBED_FILE" .qbed)

    # Define the full path for the output file
    OUTPUT_PEAK_FILE="$OUTPUT_DIR/${BASENAME}_b${BIN_SIZE}.span.peak"

    echo "----------------------------------------------------"
    echo "Starting SPAN analysis for: $BASENAME"
    echo "Output will be saved to: $OUTPUT_PEAK_FILE"
    echo "----------------------------------------------------"

    # Run the SPAN command in the background
    java $JAVA_MEM -jar "$SPAN_JAR" analyze \
        -b $BIN_SIZE \
        -kd \
        -f $FDR \
        --fragment 0 \
        -t "$QBED_FILE" \
        --format BED \
        --cs "$CHROM_SIZES" \
        -p "$OUTPUT_PEAK_FILE" &
    
    # Increment the job counter
    ((job_count++))
    
    # If the max number of jobs is reached, wait for them to finish
    if [ "$job_count" -ge "$MAX_JOBS" ]; then
        echo ""
        echo "Reached max jobs ($MAX_JOBS). Waiting for the current batch to finish..."
        wait
        job_count=0 # Reset the counter for the next batch
    fi

done

# Wait for any remaining background jobs to finish before exiting
echo ""
echo "All SPAN jobs have been started. Waiting for the final batch to complete..."
wait
echo "All jobs are finished. Results are in: $OUTPUT_DIR"
