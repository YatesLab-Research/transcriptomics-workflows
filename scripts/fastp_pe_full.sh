#!/bin/bash
set -euo pipefail

# -------- CONFIG --------
BASE_DIR="/path/to/fastq/files"
OUT_DIR="/path/to/output/pe"
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

echo "✂️ Starting fastp trimming for paired-end RNA-seq..."

for r1 in "$BASE_DIR"/*_R1*.fastq.gz "$BASE_DIR"/*_1*.fastq.gz; do
    [[ -f "$r1" ]] || continue
    r2="${r1/_R1/_R2}"
    r2="${r2/_1/_2}"

    if [[ -f "$r2" ]]; then
        sample=$(basename "$r1" | sed 's/_R1.*//; s/_1.*//')
        out_r1="$TRIMMED_DIR/${sample}_R1.trimmed.fastq.gz"
        out_r2="$TRIMMED_DIR/${sample}_R2.trimmed.fastq.gz"
        html="$LOG_DIR/${sample}_fastp.html"
        json="$LOG_DIR/${sample}_fastp.json"

        echo "🔬 Processing $sample..."
        fastp -i "$r1" -I "$r2" -o "$out_r1" -O "$out_r2" \
            --detect_adapter_for_pe \
            --qualified_quality_phred 20 \
            --length_required 50 \
            --cut_front --cut_tail \
            --trim_poly_g \
            --thread $THREADS \
            --html "$html" --json "$json"

        fastqc -t "$THREADS" -o "$FASTQC_DIR" "$out_r1" "$out_r2"
    else
        echo "⚠️  Skipping $r1 — missing R2 pair."
    fi
done

echo "📊 Running MultiQC..."
multiqc "$OUT_DIR" -o "$MULTIQC_DIR"

echo "✅ All steps complete."
echo "   Trimmed FASTQ: $TRIMMED_DIR"
echo "   FastQC reports: $FASTQC_DIR"
echo "   fastp logs: $LOG_DIR"
echo "   MultiQC summary: $MULTIQC_DIR"
