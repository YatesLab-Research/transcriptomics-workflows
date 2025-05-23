#!/bin/bash
set -euo pipefail

# -------- CONFIG --------
BAM_DIR="/path/to/hisat2_alignment_pe"  # or _se if using single-end
GTF="/path/to/annotation.gtf"
OUT_DIR="/path/to/featurecounts"
COUNT_FILE="$OUT_DIR/gene_counts.txt"
IS_PAIRED=true        # Set to false for SE reads
STRANDNESS=2          # 0: unstranded, 1: stranded, 2: reversely stranded
THREADS=8

# -------- SETUP --------
mkdir -p "$OUT_DIR"

# -------- CHECK TOOLS --------
command -v featureCounts >/dev/null || { echo "❌ featureCounts not found. Aborting."; exit 1; }

# -------- RUN FEATURECOUNTS --------
echo "📊 Running featureCounts..."

if [ "$IS_PAIRED" = true ]; then
  featureCounts -T "$THREADS" \
    -p -B -C \
    -s "$STRANDNESS" \
    -a "$GTF" -o "$COUNT_FILE" \
    "$BAM_DIR"/*.sorted.bam
else
  featureCounts -T "$THREADS" \
    -s "$STRANDNESS" \
    -a "$GTF" -o "$COUNT_FILE" \
    "$BAM_DIR"/*.sorted.bam
fi

echo "✅ featureCounts complete. Output: $COUNT_FILE"
