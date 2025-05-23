#!/bin/bash
set -euo pipefail

# Define input/output
BASE_DIR="/path/to/fastq/files"
FASTQC_OUT="/path/to/output/fastqc_pe"
MULTIQC_OUT="/path/to/output/multiqc_pe"
LOG_FILE="missing_or_mispaired_pe.log"

# Check tools
command -v fastqc >/dev/null || { echo "❌ FastQC not found"; exit 1; }
command -v multiqc >/dev/null || { echo "❌ MultiQC not found"; exit 1; }

# Prepare dirs
mkdir -p "$FASTQC_OUT" "$MULTIQC_OUT"
> "$LOG_FILE"  # empty log

echo "🚀 Starting FastQC for PE reads..."

for folder in "$BASE_DIR"/*/; do
    [ -d "$folder" ] || continue
    echo "📂 Scanning: $folder"

    for r1 in "$folder"/*_R1*.fastq.gz "$folder"/*_1*.fastq.gz; do
        [ -f "$r1" ] || continue
        base=$(basename "$r1")
        r2="${r1/_R1/_R2}"
        r2="${r2/_1/_2}"

        if [[ -f "$r2" ]]; then
            echo "🔬 Running FastQC on PE pair: $base and $(basename "$r2")"
            fastqc -o "$FASTQC_OUT" "$r1" "$r2"
        else
            echo "⚠️  Missing R2 for $base" | tee -a "$LOG_FILE"
        fi
    done
done

echo "📊 Running MultiQC summary..."
multiqc "$FASTQC_OUT" -o "$MULTIQC_OUT"

echo "✅ Completed PE QC. MultiQC report: $MULTIQC_OUT"
echo "📝 Log of missing pairs: $LOG_FILE"
