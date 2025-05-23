#!/bin/bash
set -euo pipefail

# -------- HELP --------
usage() {
  echo "Usage: $0 -c <counts_file> -m <metadata_file> -r <r_script> [-o <output_dir>]"
  echo ""
  echo "Options:"
  echo "  -c    Path to gene counts file (featureCounts output)"
  echo "  -m    Path to metadata file (TSV format)"
  echo "  -r    Path to DESeq2 R script"
  echo "  -o    Output directory (default: results/deseq2)"
  exit 1
}

# -------- PARSE ARGUMENTS --------
OUT_DIR="results/deseq2"
while getopts ":c:m:r:o:" opt; do
  case $opt in
    c) COUNTS_FILE="$OPTARG" ;;
    m) METADATA_FILE="$OPTARG" ;;
    r) R_SCRIPT="$OPTARG" ;;
    o) OUT_DIR="$OPTARG" ;;
    *) usage ;;
  esac
done

# -------- VALIDATE INPUT --------
[[ -z "${COUNTS_FILE:-}" || -z "${METADATA_FILE:-}" || -z "${R_SCRIPT:-}" ]] && usage
[[ -f "$COUNTS_FILE" ]] || { echo "❌ Counts file not found: $COUNTS_FILE"; exit 1; }
[[ -f "$METADATA_FILE" ]] || { echo "❌ Metadata file not found: $METADATA_FILE"; exit 1; }
[[ -f "$R_SCRIPT" ]] || { echo "❌ R script not found: $R_SCRIPT"; exit 1; }

# -------- SETUP --------
LOG_DIR="$OUT_DIR/logs"
mkdir -p "$OUT_DIR" "$LOG_DIR"
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
LOG_FILE="$LOG_DIR/deseq2_run_${TIMESTAMP}.log"

# -------- LOG HEADER --------
{
echo "📁 DESeq2 Pipeline Started: $(date)"
echo "🔍 Inputs:"
echo "    Counts:    $COUNTS_FILE"
echo "    Metadata:  $METADATA_FILE"
echo "    R Script:  $R_SCRIPT"
echo "    Output:    $OUT_DIR"
echo "    Log:       $LOG_FILE"

# -------- RUN R SCRIPT --------
echo "🚀 Running DESeq2 analysis..."
Rscript "$R_SCRIPT" \
  --counts "$COUNTS_FILE" \
  --metadata "$METADATA_FILE" \
  --outdir "$OUT_DIR"

echo "✅ Done: $(date)"
} 2>&1 | tee "$LOG_FILE"
