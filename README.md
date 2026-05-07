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
├── scripts/              # Analysis code
│   ├── RNA/             # RNA-seq and ATAC-seq analysis scripts
│   └── ATAC/            # ATAC-seq DiffBind analysis
├── results/             # Analysis outputs
│   ├── RNA/             # DE results, plots, tables
│   ├── ATAC/            # DA results, plots, tables
│   └── multiqc_plots/   # QC visualizations
└── metadata/            # Sample sheets

## Requirements

- R >= 4.0
- Required R packages: DESeq2, DiffBind, ChIPseeker, clusterProfiler

## Usage

Run analysis scripts in order from the project root directory.
See `scripts/RNA/run_analysis.sh` for the complete pipeline.

## Cite
