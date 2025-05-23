# RNA-seq Quality Control Pipeline (FastQC + MultiQC)

This repository contains modular scripts to run quality control (QC) on RNA-seq data using **FastQC** and **MultiQC**, with support for:

- **Paired-end (PE)** read detection
- **Single-end (SE)** read analysis
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

🧪 Usage
🧬 1. Paired-End QC
Use when your FASTQ filenames follow a paired convention (e.g.,'_R1/_R2 or _1/_2').

```bash
bash fastqc_PE.sh
```
Edit paths inside 'fastqc_PE.sh':

-'BASE_DIR' – Folder with raw FASTQ files

-'FASTQC_OUT' – Output path for FastQC reports

-'MULTIQC_OUT' – Output path for MultiQC report

What it does:

 -Detects PE files and verifies mate pairs

 -Runs FastQC on valid PE reads

 -Logs missing or unmatched files to 'missing_or_mispaired_PE.log'

 -Summarizes results via MultiQC


🔬 2. Single-End QC
Use when you have only one FASTQ file per sample.

```bash
bash fastqc_SE.sh
```
Edit paths inside 'fastqc_SE.sh':

-'BASE_DIR', 'FASTQC_OUT', 'MULTIQC_OUT', 'LOG_FILE'

What it does:

-Runs FastQC on all SE reads

-Skips '_R2' or '_2' files

-Logs ignored files to 'ignored_nonfastq_SE.log'

-Summarizes results with MultiQC

📊 Output
Each script generates:

 - FastQC Reports: HTML and zipped summaries

 - MultiQC Summary:
    -multiqc_report.html

 - Log files:

   -missing_or_mispaired_pe.log
   -ignored_nonfastq_se.log

🧼 Cleanup
To remove all outputs before a fresh run:
```bash
rm -rf /path/to/output/fastqc_pe /path/to/output/multiqc_pe
rm -rf /path/to/output/fastqc_se /path/to/output/multiqc_se
rm missing_or_mispaired_pe.log ignored_nonfastq_se.log
```

