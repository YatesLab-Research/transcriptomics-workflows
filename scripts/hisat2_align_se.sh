#!/bin/bash
set -euo pipefail

# -------- CONFIG --------
TRIMMED_DIR="/path/to/trimmed/se"
OUT_DIR="/path/to/output/hisat2_alignment_se"
HISAT2_INDEX="/path/to/hisat2/genome_index/genome"
SUMMARY_LOG="$OUT_DIR/hisat2_alignment_summary.log"
SUMMARY_CSV="$OUT_DIR/hisat2_alignment_summary.csv"
MULTIQC_LOGS="$OUT_DIR/logs"
THREADS=8

# -------- SETUP --------
mkdir -p "$OUT_DIR" "$MULTIQC_LOGS"
echo "📄 HISAT2 Alignment Summary Report" > "$SUMMARY_LOG"
echo "===================================" >> "$SUMMARY_LOG"
echo "Sample,Total,Unaligned,Aligned_once,Aligned_multi,Overall_rate" > "$SUMMARY_CSV"

# -------- CHECK TOOLS --------
for tool in hisat2 samtools; do
  command -v $tool >/dev/null || { echo "❌ $tool not found. Aborting."; exit 1; }
done

echo "🧬 Starting single-end alignment with HISAT2..."

# -------- ALIGN --------
for fq in "$TRIMMED_DIR"/*.trimmed.fastq.gz; do
    [[ -f "$fq" ]] || continue

    sample=$(basename "$fq" | sed 's/\.trimmed\.fastq\.gz//')
    bam_out="$OUT_DIR/${sample}.sorted.bam"
    log_tmp="$OUT_DIR/tmp_${sample}.log"
    log_txt="$MULTIQC_LOGS/${sample}_hisat2.log"

    echo "🔄 Aligning $sample..."
    hisat2 -x "$HISAT2_INDEX" -U "$fq" \
        --dta \
        --rna-strandness R \
        -p "$THREADS" \
        --summary-file "$log_tmp" \
        | samtools sort -@ "$THREADS" -o "$bam_out" -

    samtools index "$bam_out"

    # Append to summary log
    echo -e "\n📌 $sample" >> "$SUMMARY_LOG"
    cat "$log_tmp" >> "$SUMMARY_LOG"

    # Copy for MultiQC
    cp "$log_tmp" "$log_txt"

    # Parse to CSV
    total=$(grep "reads; of these:" "$log_tmp" | awk '{print $1}')
    unaligned=$(grep "aligned 0 times" "$log_tmp" | awk '{print $1}')
    aligned_once=$(grep "aligned exactly 1 time" "$log_tmp" | awk '{print $1}')
    aligned_multi=$(grep "aligned >1 times" "$log_tmp" | awk '{print $1}')
    rate=$(grep "overall alignment rate" "$log_tmp" | awk '{print $1}')

    echo "$sample,$total,$unaligned,$aligned_once,$aligned_multi,$rate" >> "$SUMMARY_CSV"

    rm "$log_tmp"
done

echo "✅ Single-end alignment complete."
echo "📁 BAM files: $OUT_DIR"
echo "📝 Log summary: $SUMMARY_LOG"
echo "📊 CSV summary: $SUMMARY_CSV"
echo "📦 MultiQC logs: $MULTIQC_LOGS"
