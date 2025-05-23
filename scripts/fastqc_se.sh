#!/bin/bash
set -euo pipefail

# Input/output
BASE_DIR="/path/to/fastq/files"
FASTQC_OUT="/path/to/output/fastqc_se"
MULTIQC_OUT="/path/to/output/multiqc_se"
LOG_FILE="ignored_nonfastq_se.log"

# Check tools
command -v fastqc >/dev/null || { echo "❌ FastQC not found"; exit 1; }
command -v multiqc >/dev/null || { echo "❌ MultiQC not found"; exit 1; }

mkdir -p "$FASTQC_OUT" "$MULTIQC_OUT"
> "$LOG_FILE"

echo "🚀 Starting FastQC for SE reads..."

for folder in "$BASE_DIR"/*/; do
    [ -d "$folder" ] || continue
    echo "📂 Scanning: $folder"

    for file in "$folder"/*.fastq.gz; do
        [[ -f "$file" ]] || continue

        filename=$(basename "$file")

        # Skip R2 or known PE files
        if [[ "$filename" == *_R2* || "$filename" == *_2* ]]; then
            echo "🔁 Skipping paired-end file: $filename"
            continue
        fi

        # Check for valid FASTQ extension
        if [[ "$filename" =~ \.fastq\.gz$ ]]; then
            echo "🔬 Running FastQC on: $filename"
            fastqc -o "$FASTQC_OUT" "$file"
        else
            echo "⚠️  Ignored non-fastq file: $filename" | tee -a "$LOG_FILE"
        fi
    done
done

echo "📊 Running MultiQC summary..."
multiqc "$FASTQC_OUT" -o "$MULTIQC_OUT"

echo "✅ Completed SE QC. MultiQC report: $MULTIQC_OUT"
echo "📝 Log of ignored files: $LOG_FILE"
