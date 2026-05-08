# RNA-seq and ATAC-seq Analysis Pipeline

Analysis code for integrated RNA-seq and ATAC-seq differential analysis.

## Data Availability
##
Raw sequencing data will be deposited to NCBI GEO.

## Overview

This repository contains analysis scripts for:
- RNA-seq differential expression (DESeq2)
- ATAC-seq differential accessibility (DiffBind)  
- Integrated multi-omics analysis

## Directory Structure
<<<<<<< HEAD
=======

```text
├── scripts/              # Analysis code 
│   |── RNA/             # RNA-seq and ATAC-seq analysis scripts 
│   └── ATAC/            # ATAC-seq DiffBind analysis 
├── results/             # Analysis outputs 
│   ├── RNA/             # DE results, plots, tables 
│   ├── ATAC/            # DA results, plots, tables 
│   └── multiqc_plots/   # QC visualizations 
└── metadata/            # Sample sheets 
```
>>>>>>> e3684a7fe48dbd2ac8d785a328db89f6ad8d44b2

```text
├── scripts/              # Analysis code 
│   |── RNA/             # RNA-seq and ATAC-seq analysis scripts 
│   └── ATAC/            # ATAC-seq DiffBind analysis 
├── results/             # Analysis outputs 
│   ├── RNA/             # DE results, plots, tables 
│   ├── ATAC/            # DA results, plots, tables 
│   └── multiqc_plots/   # QC visualizations 
└── metadata/            # Sample sheets 
```

## Software Versions

### R Environment
- R: 4.2.3

### Differential Expression Analysis
- DESeq2: 1.38.3
- edgeR: 3.40.2
- limma: 3.54.2
- RUVSeq: 1.32.0

### Differential Accessibility Analysis
- DiffBind: 3.20.0
- edgeR: 4.8.2
- limma: 3.66.0

### Peak Annotation
- ChIPseeker: 1.46.1

### Functional Enrichment Analysis
- clusterProfiler: 4.6.2
- enrichplot: 1.18.4
- org.Hs.eg.db: 3.16.0
- GO.db: 3.16.0


### Data Visualization
- ggplot2: 4.0.2
- dplyr: 1.2.0
- tidyr: 1.3.2
- pheatmap: 1.0.13
- RColorBrewer: 1.1-3

## Usage

Run analysis scripts in order from the project root directory.
See `scripts/RNA/run_analysis.sh` for the complete pipeline.

## Cite


