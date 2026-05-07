#!/usr/bin/env Rscript
# Package installation for RNA-seq and ATAC-seq analysis

# CRAN packages
cran_pkgs <- c(
  "tidyverse", "janitor", "ggrepel", "ggpubr",
  "pheatmap", "patchwork", "RColorBrewer"
)

for (pkg in cran_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}

# BiocManager
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", repos = "https://cloud.r-project.org")

# Bioconductor packages
bioc_pkgs <- c(
  "DESeq2", "RUVSeq", "EDASeq", "limma", "ggplot2",
  "clusterProfiler", "org.Hs.eg.db", "enrichplot",
  "DiffBind", "ChIPseeker", "TxDb.Hsapiens.UCSC.hg38.knownGene",
  "GenomicRanges", "rtracklayer", "apeglm", "ashr"
)

for (pkg in bioc_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
  }
}

sessionInfo()
