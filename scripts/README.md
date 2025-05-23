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
Edit paths inside 'fastqc_pe.sh':

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
