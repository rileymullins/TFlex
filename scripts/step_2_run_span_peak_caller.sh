#!/bin/bash
# step_2_run_span_peak_caller.sh

# ==============================================================================
# Default Configuration
# ==============================================================================
# 
# --- SPAN Parameters for ---
JAVA_MEM="-Xmx16G"
BIN_SIZE=50
FDR="0.01"

# ==============================================================================
# Command-Line Argument Parsing
# ==============================================================================

# Function to display script usage
usage() {
    echo "Usage: $0 [-j SPAN_JAR] [-i INPUT_DIR] [-o OUTPUT_DIR] [-c CHROM_SIZES] [-b BIN_SIZE] [-f FDR]"
    echo " "
    echo "Options:"
    echo "  -j      Path to the SPAN jar file."
    echo "  -i      Directory containing input .qbed files."
    echo "  -o      Directory for output peak files."
    echo "  -c      Path to the chromosome sizes file."
    echo "  -b      Bin size for SPAN analysis."
    echo "  -f      False Discovery Rate (FDR) cutoff. (Default: $FDR)"
    echo "  -h      Display this help message."
    exit 1
}

# Parse the options
while getopts "j:i:o:c:b:f:h" opt; do
    case ${opt} in
        j ) SPAN_JAR="$OPTARG" ;;
        i ) INPUT_DIR="$OPTARG" ;;
        o ) OUTPUT_DIR="$OPTARG" ;;
        c ) CHROM_SIZES="$OPTARG" ;;
        b ) BIN_SIZE="$OPTARG" ;;
        f ) FDR="$OPTARG" ;;
        h ) usage ;;
        \? ) usage ;;
    esac
done

# ==============================================================================
# Script Logic
# ==============================================================================

# Create the output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Find all .qbed files in the input directory and loop through them
for QBED_FILE in "$INPUT_DIR"/*.qbed; do
    # Check if files exist to avoid errors with empty directories
    [ -e "$QBED_FILE" ] || continue

    # Get the base filename without the path and extension
    BASENAME=$(basename "$QBED_FILE" .qbed)

    # Define the full path for the output file
    OUTPUT_PEAK_FILE="$OUTPUT_DIR/${BASENAME}_b${BIN_SIZE}.span.peak"

    echo "----------------------------------------------------"
    echo "Starting SPAN analysis for: $BASENAME"
    echo "Output will be saved to: $OUTPUT_PEAK_FILE"
    echo "----------------------------------------------------"

    # Run the SPAN command sequentially. Recommended bin size and FDR value are given.
    java $JAVA_MEM -jar "$SPAN_JAR" analyze \
        -b 50 \
        -kd \
        -f $FDR \
        --fragment 0 \
        -t "$QBED_FILE" \
        --format BED \
        --cs "$CHROM_SIZES" \
        -p "$OUTPUT_PEAK_FILE"
done

echo ""
echo "All jobs are finished. Results are in: $OUTPUT_DIR"
