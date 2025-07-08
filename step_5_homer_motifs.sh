#!/bin/bash

# === HOMER Motif Discovery Batch Script ===
#
# Runs HOMER findMotifsGenome.pl on each BED file produced by your R pipeline.
# Each BED file represents peaks for a TF group (e.g. JUN_significant_peaks.bed).
#
# === Usage ===
#   bash run_homer_motifs.sh <BED_DIR> <GENOME> [<NUM_WORKERS>]
#
# Arguments:
#   <BED_DIR>     : Directory containing *_significant_peaks.bed files
#   <GENOME>      : Genome assembly (e.g. hg38, mm10)
#   [<NUM_WORKERS>]: (Optional) Number of parallel jobs to run.
#                    If omitted, defaults to half your CPU cores (min 1).
#
# Example:
#   bash run_homer_motifs.sh /path/to/MotifAnalysis hg38 4
#
# === Requirements ===
# - HOMER (findMotifsGenome.pl must be in PATH)
# - GNU parallel (optional; script will use if available)
#
# === Output ===
# - Motif results saved to: <BED_DIR>/HOMER_Motif_Output/<TF>_motif_output/
#

set -euo pipefail

# === Parse arguments ===
if [[ $# -lt 2 || $# -gt 3 ]]; then
  echo "Usage: $0 <BED_DIR> <GENOME> [<NUM_WORKERS>]"
  echo "Example: $0 /path/to/MotifAnalysis hg38 4"
  exit 1
fi

BED_DIR="$1"
GENOME="$2"
OUT_DIR="${BED_DIR}/HOMER_Motif_Output"

# If provided, use user-defined workers; else leave blank to trigger auto mode
NUM_WORKERS="${3:-}"

mkdir -p "$OUT_DIR"

# === Check HOMER ===
if ! command -v findMotifsGenome.pl &>/dev/null; then
  echo "Error: findMotifsGenome.pl not found in PATH. Ensure HOMER is installed and in your PATH." >&2
  exit 1
fi

# === Function to run HOMER ===
run_homer() {
  local bed_file="$1"
  local tf_name
  tf_name=$(basename "$bed_file" | sed 's/_significant_peaks\.bed$//')
  local tf_outdir="${OUT_DIR}/${tf_name}_motif_output"
  findMotifsGenome.pl "$bed_file" "$GENOME" "$tf_outdir" -size given -mask
}

export -f run_homer
export OUT_DIR
export GENOME

# === Collect BED files ===
bed_files=("${BED_DIR}"/*_significant_peaks.bed)
num_files=${#bed_files[@]}

if [[ $num_files -eq 0 ]]; then
  echo "No *_significant_peaks.bed files found in $BED_DIR"
  exit 1
fi

echo "Found $num_files BED files. Starting HOMER motif analysis..."

# === Run HOMER ===
if command -v parallel &>/dev/null; then
  # Determine number of jobs
  if [[ -n "$NUM_WORKERS" ]]; then
    jobs="$NUM_WORKERS"
  else
    cores=$(getconf _NPROCESSORS_ONLN)
    jobs=$(( cores / 2 ))
    (( jobs < 1 )) && jobs=1
  fi

  echo "Using $jobs parallel jobs for HOMER runs."
  printf "%s\n" "${bed_files[@]}" | parallel --bar --jobs "$jobs" run_homer {}
else
  echo "GNU parallel not found, running background jobs with manual progress..."
  completed=0
  total=$num_files
  pids=()

  show_progress() {
    local current=$1
    local total=$2
    local width=40
    local percent=$(( 100 * current / total ))
    local filled=$(( width * current / total ))
    local empty=$(( width - filled ))
    printf "\r["
    for ((i=0; i<filled; i++)); do printf "#"; done
    for ((i=0; i<empty; i++)); do printf " "; done
    printf "] %d/%d (%d%%)" "$current" "$total" "$percent"
  }

  for bed_file in "${bed_files[@]}"; do
    run_homer "$bed_file" &
    pids+=($!)
  done

  for pid in "${pids[@]}"; do
    wait "$pid"
    ((completed++))
    show_progress $completed $total
  done
  echo
fi

echo "All HOMER motif analyses complete."
echo "Results stored in: $OUT_DIR"
