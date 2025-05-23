# RNA-seq Data Analysis Workflow 
---

# RNA-seq Quality Control (FastQC + MultiQC)

This repository contains modular scripts (```fastqc_pe.sh ``` & ```fastqc_se.sh ```) to run quality control (QC) on RNA-seq data using **FastQC** and **MultiQC**, with support for:

- **Paired-end (pe)** read detection
- **Single-end (se)** read analysis
- **Automated logging** of missing or mispaired files

---

## 🔧 Requirements

Install and ensure the following tools are in your `$PATH`:

- [`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- [`MultiQC`](https://multiqc.info/)

Install via Conda (recommended):

```bash
conda install -c bioconda fastqc multiqc
```

## 🧪 **Usage**

## 🧬 1. Paired-End QC
Use when your FASTQ filenames follow a paired convention (e.g.,```_R1/_R2``` or```_1/_2```).

```bash
bash fastqc_pe.sh
```
**Edit paths inside 'fastqc_pe.sh':**

 -```BASE_DIR``` – Folder with raw FASTQ files
 
 -```FASTQC_OUT``` – Output path for FastQC reports
 
 -```MULTIQC_OUT``` – Output path for MultiQC report

**What it does:**

  -Detects PE files and verifies mate pairs
  
  -Runs FastQC on valid PE reads
  
  -Logs missing or unmatched files to ```missing_or_mispaired_pe.log```
  
  -Summarizes results via MultiQC

## 🔬 2. Single-End QC
Use when you have only one FASTQ file per sample.

```bash
bash fastqc_se.sh
```
Edit paths inside 'fastqc_se.sh':

-```BASE_DIR```, ```FASTQC_OUT```, ```MULTIQC_OUT```, ```LOG_FILE```

**What it does:**
-Runs FastQC on all SE reads
-Skips ```_R2``` or ```_2``` files
-Logs ignored files to ```ignored_nonfastq_se.log```
-Summarizes results with MultiQC

## 📊  Output
Each script generates:

 - FastQC Reports: HTML and zipped summaries
 - MultiQC Summary: ```multiqc_report.html```
 - Log files:
   -```missing_or_mispaired_pe.log```
   -```ignored_nonfastq_se.log```

## 🧼 Cleanup
To remove all outputs before a fresh run:
```bash
rm -rf /path/to/output/fastqc_se /path/to/output/multiqc_se
rm -rf /path/to/output/fastqc_se /path/to/output/multiqc_se
rm missing_or_mispaired_se.log ignored_nonfastq_se.log
```
> The RNA-seq data quality control was performed using custom Bash scripts built around FastQC and MultiQC.

## 🔎 Interpreting FastQC to Guide Trimming

Before trimming with `fastp` or any other tool, it's good practice to run **initial FastQC** on raw reads to assess:

| FastQC Flag                   | How to Respond in `fastp`                              |
|------------------------------|---------------------------------------------------------|
| Adapter Content (❌)         | Use `--detect_adapter_for_pe` or specify adapters      |
| Per Base Quality Drop (❌)   | Enable `--cut_front` and `--cut_tail` trimming         |
| Poly-G Tail (NovaSeq data)   | Add `--trim_poly_g`                                     |
| Short Reads (Post-QC)        | Use `--length_required 50` to discard <50bp fragments  |
| Overrepresented Sequences    | Consider removing with custom adapter lists or filtering |

---

# ✂️ RNA-seq Trimming and QC with fastp
This repository includes two scripts (`fastp_pe_full.sh` & `fastp_se_full.sh`) for trimming RNA-seq reads with fastp, followed by quality control with FastQC, and summarization with MultiQC.

## 🔧 Requirements
Install tools (e.g. via Conda):
```bash
conda install -c bioconda fastp fastqc multiqc
```

> After identifying issues with **FastQC** from the first step above, adjust the `fastp` settings accordingly (these are built into the provided scripts).

## 🧬 1. Paired-End Reads
Script: ```fastp_pe_full.sh```

```bash
bash fastp_pe_full.sh
```
**What it does:**

-Detects paired-end reads (```_R1```/```_R2``` or ```_1```/```_2```)
-Trims adapters and low-quality ends
-Filters short reads (length < 50 bp)
-Runs FastQC on trimmed output
-Generates a MultiQC report

**Edit the script to set:**
-```BASE_DIR```: Input FASTQ folder
-```OUT_DIR```: Output base folder

## 🔬 2. Single-End Reads
Script: fastp_se_full.sh

```bash
bash fastp_se_full.sh
```
**What it does:**

-Detects all non-R2 ```.fastq.gz``` files
-Trims adapters and poly-G tails
-Filters low-quality and short reads
-Runs FastQC and MultiQC as above

📁 Output Structure

```graphql
output_dir/
├── trimmed/         # trimmed FASTQ files
├── fastqc/          # FastQC reports for trimmed files
├── fastp_logs/      # fastp HTML and JSON reports
└── multiqc/         # combined MultiQC report
```

## ✅ Example Use

```bash
# Paired-end trimming
./fastp_pe_full.sh

# Single-end trimming
./fastp_se_full.sh
```
---

# 🧬 RNA-seq Alignment with HISAT2

This section describes how to align trimmed RNA-seq reads to a reference genome using HISAT2. Both paired-end and single-end pipelines are supported, with alignment statistics logged for downstream analysis and MultiQC integration.

## Requirements

Install tools via Conda:

```bash
conda install -c bioconda hisat2 samtools multiqc
```

## Prepare the Reference Genome

Download a reference genome FASTA and annotation (GTF) file, e.g., from GENCODE.
Example:

```bash
# Download
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_45/GRCh38.primary_assembly.genome.fa.gz
gunzip GRCh38.primary_assembly.genome.fa.gz

# Build HISAT2 index
hisat2-build -p 8 GRCh38.primary_assembly.genome.fa genome_index/genome
```
This will create a ```genome_index/``` folder containing ```.ht2``` files. Use the path to the index base (```genome_index/genome```) in the scripts below.


## 🔁 1. Paired-End Alignment

Script: ```hisat2_align_pe.sh```

```bash
bash hisat2_align_pe.sh
```

**What it does:**
-Aligns PE reads (```_R1```/```_R2```) from fastp
-Outputs sorted ```.bam``` and ```.bai``` files
-Logs alignment statistics for each sample
-Compiles a combined ```.log``` and ```.csv``` file
-Supports MultiQC log detection

**Edit inside script:**
-```TRIMMED_DIR```: path to your PE trimmed FASTQs
-```HISAT2_INDEX```: base path to HISAT2 index
-```OUT_DIR```: where outputs go

## 🔁 2. Single-End Alignment
Script: hisat2_align_se.sh

```bash
bash hisat2_align_se.sh
```

**Edit inside script:**

-```TRIMMED_DIR```, ```HISAT2_INDEX```, and ```OUT_DIR```

## 📊 Output Structure

Each alignment script will generate:

```
hisat2_alignment_*/                     # PE or SE
├── sample1.sorted.bam                  # Aligned, sorted reads
├── sample1.sorted.bam.bai              # BAM index
├── hisat2_alignment_summary.csv        # CSV for downstream plotting
├── hisat2_alignment_summary.log        # Combined readable summary
├── logs/                               # MultiQC-compatible logs
│   ├── sample1_hisat2.log
│   └── sample2_hisat2.log
```

## 📈 MultiQC Integration

After alignment, run:

```bash
multiqc /path/to/hisat2_alignment_pe  -o /path/to/hisat2_alignment_pe/multiqc
multiqc /path/to/hisat2_alignment_se  -o /path/to/hisat2_alignment_se/multiqc
```

## Bonus

## 🧬 HISAT2 RNA-seq Alignment: Key Parameters

| Parameter             | Description                                   | Paired-End | Single-End | Recommended Value / Notes                          |
|-----------------------|-----------------------------------------------|------------|------------|----------------------------------------------------|
| `-x`                  | Path to HISAT2 index basename                 | ✅         | ✅         | e.g. `genome_index/genome` (omit `.ht2` extension) |
| `-1/-2`               | Input R1 and R2 reads                         | ✅         | ❌         | `/path/sample_R1.fastq.gz`, `/path/sample_R2.fastq.gz` |
| `-U`                  | Input for single-end reads                    | ❌         | ✅         | `/path/sample.fastq.gz`                            |
| `-p`                  | Number of threads                             | ✅         | ✅         | 4–16 (depending on machine)                        |
| `--dta`               | Optimize for transcript assembly              | ✅         | ✅         | Always recommended for RNA-seq                     |
| `--rna-strandness`    | Strand-specific protocol setting              | ✅         | ✅         | `RF`, `FR`, `R`, `F` (depends on library prep)     |
| `--summary-file`      | Write per-sample alignment stats              | ✅         | ✅         | Used for logging & MultiQC                         |
| `| samtools sort`     | Pipe to sort BAM output                       | ✅         | ✅         | Required for downstream tools                      |
| `samtools index`      | Create BAM index (`.bai`)                     | ✅         | ✅         | Required for IGV and counting                      |


## 📌 Common Values for `--rna-strandness`

| Library Type (Protocol)   | `--rna-strandness` Value | Notes                              |
|---------------------------|--------------------------|------------------------------------|
| TruSeq Stranded / dUTP    | `RF`                     | Most common in Illumina RNA-seq    |
| Illumina Unstranded       | *omit flag*              | Default                            |
| Reverse (3′ → 5′)         | `R`                      | Use with single-end libraries      |
| Forward (5′ → 3′)         | `F`                      | Use with single-end libraries      |

## 🧪 Output Files to Expect

| Output Type       | Description                                         |
|-------------------|-----------------------------------------------------|
| `sample.sorted.bam` | Sorted alignments ready for quantification         |
| `sample.sorted.bam.bai` | BAM index for IGV and featureCounts             |
| `sample_hisat2.log` | Per-sample HISAT2 log (alignment stats)           |
| `hisat2_alignment_summary.csv` | Combined summary table for downstream plotting |
| `hisat2_alignment_summary.log` | Human-readable combined summary          |

---
# 🧬 RNA-seq Gene Quantification with featureCounts

This module performs gene-level quantification from aligned RNA-seq reads (BAM files) using [`featureCounts`](http://bioinf.wehi.edu.au/featureCounts/), part of the **Subread** package.

---

## 📦 Requirements

Install via Conda:

```bash
conda install -c bioconda subread
```
## 🗂 Input Requirements
-Sorted BAM files from HISAT2 (PE or SE)
-GTF annotation file matching the reference genome

> 🔗 Download GTF from GENCODE: https://www.gencodegenes.org

Example:

```bash
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_45/gencode.v45.annotation.gtf.gz
gunzip gencode.v45.annotation.gtf.gz
```

## 🚀 How to Run
Edit and run the provided script:

```bash
bash run_featurecounts.sh
```

##📄 Output
-`gene_counts.txt`: Count matrix for downstream tools (e.g. DESeq2)
-Columns:
  -Gene ID
  -Read counts per sample

## 📌 Best Practices
-Always match the GTF file to the reference genome used during alignment.
-Use the correct strandness (`-s`) setting (	`0`, `1`, or `2`) — check with your sequencing provider or tools like `RSeQC`.
-Use sorted and indexed BAMs from HISAT2.

---

# 🧭 What's Next: Downstream Analysis in R

This marks the **end of the Bash-based RNA-seq processing pipeline**.

All upstream steps — quality control, trimming, alignment, and gene-level quantification — have been completed.

The **next stage** involves downstream statistical and biological analysis in **R**.

---

## 📊 Downstream Analysis in R (Suggestions)

| Goal                        | R Tool / Package     | Description                                               |
|-----------------------------|----------------------|-----------------------------------------------------------|
| Differential expression     | **DESeq2** / **edgeR** | Identify up/down-regulated genes between sample groups   |
| Visualization               | **ggplot2**, **pheatmap**, **ComplexHeatmap** | Volcano plots, heatmaps, PCA                              |
| Data normalization          | **DESeq2**, **limma-voom** | Normalize raw counts to account for library size          |
| Gene annotation             | **biomaRt**, **org.Hs.eg.db** | Map Ensembl IDs to gene symbols, pathways, etc.         |
| Pathway enrichment          | **clusterProfiler**, **fgsea**, **GSEA** | Explore functional implications of DE genes              |
| Quality assessment (PCA, clustering) | **DESeq2**, **Rtsne**, **PCAtools** | Detect batch effects, sample relationships               |

---

# 🎯 DESeq2 Differential Expression Analysis Pipeline

This pipeline performs differential gene expression analysis using **DESeq2** with support for:
- 🧬 Batch effect correction
- 🔁 Custom contrast group comparisons
- 🏷 Gene annotation via biomaRt
- 📈 Publication-quality plots (MA, Volcano, Heatmap, PCA)
- 🧮 Summary of DE gene counts

---

## ⚙️ Requirements

Install required packages (first-time only):

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("DESeq2", "EnhancedVolcano", "biomaRt", "pheatmap", "RColorBrewer", "ggplot2"))
```

## 📂 File Structure

```bash
project/
├── data/
│   ├── gene_counts.txt       # from featureCounts
│   └── metadata.tsv          # sample metadata (with 'condition' and 'batch')
├── deseq2_analysis_advanced.R
├── run_deseq2.sh
└── results/
    └── <your_run_name>/
        ├── deseq2_results_treatment_vs_control.csv
        ├── deseq2_results_annotated_treatment_vs_control.csv
        ├── de_gene_summary.txt
        └── figures/
            ├── ma_plot.png/.pdf
            ├── volcano_plot.png/.pdf
            ├── heatmap_top30.png/.pdf
            └── pca_plot.png/.pdf
```

## 🚀 How to Run

```bash
bash run_deseq2.sh \
  -c data/gene_counts.txt \
  -m data/metadata.tsv \
  -r deseq2_analysis_advanced.R \
  -o results/my_analysis_run \
  --groupA treatment \
  --groupB control
```

###❗ Important Reminders
> The values you use for `--groupA` and `--groupB` must match exactly the `condition` values in your `metadata.tsv` file (case-sensitive).

## 📊 Output Summary

| File                                      | Description                                                    |
| ----------------------------------------- | -------------------------------------------------------------- |
| `deseq2_results_<A>_vs_<B>.csv`           | DESeq2 base results                                            |
| `deseq2_results_annotated_<A>_vs_<B>.csv` | With HGNC symbols & gene descriptions                          |
| `de_gene_summary.txt`                     | Count of upregulated, downregulated, and not significant genes |
| `figures/*.png/.pdf`                      | MA, Volcano, Heatmap, PCA — publication-ready                  |


## 📝 Metadata Format

Tab-separated .tsv with headers:

```tsv
sample	condition	batch
Sample1	control	Batch1
Sample2	control	Batch1
Sample3	control	Batch2
Sample4	treatment	Batch1
Sample5	treatment	Batch2
Sample6	treatment	Batch2
```

## 📌 Notes
-Adjust ```--groupA``` and ```--groupB``` to test other conditions such as `Healthy` vs `disease` etc.

-Supports batch correction automatically via ```design = ~ batch + condition```.

-Set input/output paths flexibly via CLI.

-All visual outputs saved in both PNG and PDF formats.

