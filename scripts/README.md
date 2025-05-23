# RNA-seq Quality Control Pipeline (FastQC + MultiQC)

This repository contains modular scripts to run quality control (QC) on RNA-seq data using **FastQC** and **MultiQC**, with support for:

- **Paired-end (PE)** read detection
- **Single-end (SE)** read analysis
- **Automated logging** of missing or mispaired files

---

## 📁 Folder Structure

.
├── fastqc_pe.sh # PE QC script with auto-detection and logging
├── fastqc_se.sh # SE QC script with basic filtering and logging
├── missing_or_mispaired_pe.log # Log of PE reads missing their mate
├── ignored_nonfastq_se.log # Log of SE files that were skipped
└── README.md


---

## 🔧 Requirements

Install and ensure the following tools are in your `$PATH`:

- [`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- [`MultiQC`](https://multiqc.info/)

Install via Conda (recommended):

```bash
conda install -c bioconda fastqc multiqc
```
Edit paths inside fastqc_se.sh:

-BASE_DIR, FASTQC_OUT, MULTIQC_OUT, LOG_FILE

What it does:

-Runs FastQC on all SE reads

-Skips _R2 or _2 files

-Logs ignored files to ignored_nonfastq_se.log

-Summarizes results with MultiQC

📊 Output
Each script generates:

 - FastQC Reports: HTML and zipped summaries

 - MultiQC Summary:
 -      -multiqc_report.html

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

