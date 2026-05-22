# RNA-seq and ATAC-seq Analysis Pipeline

Analysis code for the Kabuki syndrome type I fibroblast multi-omics dataset.

## Data Availability

Raw sequencing data and processed files are deposited in NCBI GEO:
- ATAC-seq: [GSE330760](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE330760)
- RNA-seq: [GSE330762](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE330762)

Data are under embargo and will be publicly released upon publication.

## Overview

This repository contains analysis scripts for:
- RNA-seq differential expression (DESeq2 + RUVSeq)
- ATAC-seq differential accessibility (DiffBind + DESeq2 + RUVSeq)
- Integrated RNA–ATAC multi-omics analysis

## Directory Structure
├── scripts/
│   ├── install_packages.R     # 0. Install dependencies
│   ├── run_rnaseq.sh          # 1. nf-core/rnaseq pipeline
│   ├── run_atacseq.sh         # 2. nf-core/atacseq pipeline
│   ├── run_diffbind.sh        # 3. DiffBind peak counting
│   ├── run_subset.sh          # 4. BAM subsetting for IGV
│   ├── 00_config.R            # Shared parameters and paths
│   ├── B01_QC_collapse.R      # QC and replicate collapsing
│   ├── B02_RNA_DE.R           # Differential expression
│   ├── B03_ATAC_DA.R          # Differential accessibility
│   ├── B04_integration.R      # RNA–ATAC integration
│   ├── diffbind_analysis.R    # DiffBind consensus peaks
│   ├── sample_sheet.csv       # Sample metadata
│   └── targets.bed            # Target regions for BAM subsetting
├── metadata/
│   ├── samplesheet_RNA.csv    # nf-core/rnaseq input
│   └── samplesheet_ATAC.csv   # nf-core/atacseq input
└── README.md

## Usage

Run scripts in order from the project root directory:

1. `Rscript scripts/install_packages.R`
2. `bash scripts/run_rnaseq.sh`
3. `bash scripts/run_atacseq.sh`
4. `bash scripts/run_diffbind.sh`
5. `Rscript scripts/B01_QC_collapse.R`
6. `Rscript scripts/B02_RNA_DE.R`
7. `Rscript scripts/B03_ATAC_DA.R`
8. `Rscript scripts/B04_integration.R`

Figure assembly scripts are not included; all source data required to reproduce figures are available via GEO.

## Software Versions

| Tool | Version |
|------|---------|
| R | 4.2.3 |
| DESeq2 | 1.38.3 |
| RUVSeq | 1.32.0 |
| limma | 3.54.2 |
| DiffBind | 3.20.0 |
| ChIPseeker | 1.34.1 |
| clusterProfiler | 4.6.2 |
| nf-core/rnaseq | 3.14.0 |
| nf-core/atacseq | 2.1.2 |

## Citation

Danvin E. et al. Transcriptomic and chromatin accessibility dataset from Kabuki syndrome type I skin fibroblasts. *Scientific Data* (2026).
