#!/bin/bash
set -euo pipefail

# -------- CONFIG --------
BASE_DIR="/path/to/fastq/files"
OUT_DIR="/path/to/output/se"
TRIMMED_DIR="$OUT_DIR/trimmed"
FASTQC_DIR="$OUT_DIR/fastqc"
LOG_DIR="$OUT_DIR/fastp_logs"
MULTIQC_DIR="$OUT_DIR/multiqc"
THREADS=8

# -------- PREP --------
mkdir -p "$TRIMMED_DIR" "$FASTQC_DIR" "$LOG_DIR" "$MULTIQC_DIR"

# -------- CHECK TOOLS --------
for tool in fastp fastqc multiqc; do
  command -v $tool >/dev/null || { echo "❌ $tool not found in PATH. Aborting."; exit 1; }
done

echo "✂️ Starting fastp trimming for single-end RNA-seq..."

for f in "$BASE_DIR"/*.fastq.gz; do
    [[ -f "$f" ]] || continue
    fname=$(basename "$f")

    # Skip if it's part of a PE pair
    [[ "$fname" =~ _R2|_2 ]] && continue

    sample=$(basename "$f" | sed 's/.fastq.gz//')
    out_trim="$TRIMMED_DIR/${sample}.trimmed.fastq.gz"
    html="$LOG_DIR/${sample}_fastp.html"
    json="$LOG_DIR/${sample}_fastp.json"

    echo "🔬 Processing $sample..."
    fastp -i "$f" -o "$out_trim" \
          --qualified_quality_phred 20 \
          --length_required 50 \
          --cut_front --cut_tail \
          --trim_poly_g \
          --thread $THREADS \
          --html "$html" --json "$json"

    fastqc -t "$THREADS" -o "$FASTQC_DIR" "$out_trim"
done

echo "📊 Running MultiQC..."
multiqc "$OUT_DIR" -o "$MULTIQC_DIR"

echo "✅ All steps complete."
echo "   Trimmed FASTQ: $TRIMMED_DIR"
echo "   FastQC reports: $FASTQC_DIR"
echo "   fastp logs: $LOG_DIR"
echo "   MultiQC summary: $MULTIQC_DIR"
