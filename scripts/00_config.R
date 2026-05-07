# Configuration for RNA-seq and ATAC-seq analysis

suppressPackageStartupMessages({
  library(tidyverse)
  library(janitor)
})

PROJECT <- "."

# File paths
PATHS <- list(
  sample_sheet     = file.path(PROJECT, "scripts/sample_sheet.csv"),
  counts_raw       = file.path(PROJECT, "pipeline/RNA/results/star_salmon/salmon.merged.gene_counts.tsv"),
  atac_diffbind    = file.path(PROJECT, "results/ATAC/DiffBind/Background_norm/diffbind_analyzed.rds"),
  counts_collapsed = file.path(PROJECT, "results/RNA/results/tables/B01_collapsed_counts.tsv"),
  meta_collapsed   = file.path(PROJECT, "results/RNA/results/tables/B01_collapsed_metadata.csv"),
  gene_map         = file.path(PROJECT, "results/RNA/results/tables/B01_gene_map.tsv"),
  rna_de           = file.path(PROJECT, "results/RNA/results/tables/B02_RNA_DE_results.csv"),
  rna_go           = file.path(PROJECT, "results/RNA/results/tables/B02_RNA_GO_BP.csv"),
  dds_rna          = file.path(PROJECT, "results/RNA/results/rds/B02_dds_RNA_final.rds"),
  atac_da          = file.path(PROJECT, "results/RNA/results/tables/B03_ATAC_DA_peaks.csv"),
  atac_annotated   = file.path(PROJECT, "results/RNA/results/tables/B03_ATAC_DA_peaks_annotated.csv"),
  dds_atac         = file.path(PROJECT, "results/RNA/results/rds/B03_dds_ATAC_final.rds"),
  integration      = file.path(PROJECT, "results/RNA/results/tables/B04_concordant_genes.csv"),
  plots            = file.path(PROJECT, "results/RNA/results/plots"),
  tables           = file.path(PROJECT, "results/RNA/results/tables"),
  rds              = file.path(PROJECT, "results/RNA/results/rds"),
  session          = file.path(PROJECT, "results/RNA/results/session_info")
)

# Load cohort metadata filtered by assay type
load_cohort <- function(assay = c("RNA", "ATAC")) {
  assay <- match.arg(assay)
  ss <- read.csv(PATHS$sample_sheet, stringsAsFactors = FALSE) %>%
    janitor::clean_names()
  assay_string <- if (assay == "RNA") "RNA-seq" else "ATAC"
  ss <- ss %>%
    dplyr::filter(grepl(assay_string, assay_performed, fixed = TRUE)) %>%
    dplyr::mutate(
      status = factor(status, levels = c("control", "KS_I")),
      sex    = factor(sex),
      batch  = factor(batch)
    ) %>%
    dplyr::distinct(gdb, .keep_all = TRUE) %>%
    dplyr::arrange(gdb)
  rownames(ss) <- ss$gdb
  ss
}

# Analysis parameters
K_FACTOR        <- 2          # RUVSeq k factors
N_EMPIRICAL     <- 5000       # Number of empirical genes for RUVg
P_CUTOFF        <- 0.05       # FDR threshold
LFC_CUTOFF      <- 1.0        # Strict log2 fold change threshold
LFC_CUTOFF_SOFT <- 0.585      # Relaxed log2 fold change threshold
GO_ONT          <- "BP"       # Gene Ontology: Biological Process
GO_SIMPLIFY     <- 0.7        # GO term simplification cutoff

# Color schemes
COLORS_STATUS  <- c(control = "grey80", KS_I = "firebrick")
COLORS_SEX     <- c(M = "#56B4E9", F = "#E69F00")
COLORS_BATCH   <- c("84" = "#1B7837", "185" = "#762A83")
COLORS_VOLCANO <- c(Upregulated = "firebrick3", Downregulated = "navy",
                    `Not Significant` = "grey85")

# Initialize output directories
init_dirs <- function() {
  for (d in c(PATHS$plots, PATHS$tables, PATHS$rds, PATHS$session))
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# Save ggplot object to PDF
save_plot <- function(plot_obj, name, width = 8, height = 8) {
  init_dirs()
  path <- file.path(PATHS$plots, paste0(name, ".pdf"))
  ggplot2::ggsave(path, plot_obj, device = "pdf", width = width, height = height)
  invisible(plot_obj)
}

# Save pheatmap directly to PDF
save_heatmap <- function(mat, anno_col, anno_colors, title, name,
                         width = 8, height = 10,
                         show_rownames = TRUE, fontsize_row = 6, ...) {
  init_dirs()
  path <- file.path(PATHS$plots, paste0(name, ".pdf"))
  pheatmap::pheatmap(
    mat,
    annotation_col    = anno_col,
    annotation_colors = anno_colors,
    color             = colorRampPalette(c("navy", "white", "firebrick"))(100),
    scale             = "row",
    main              = title,
    filename          = path,
    width             = width,
    height            = height,
    show_rownames     = show_rownames,
    fontsize_row      = fontsize_row,
    ...
  )
}

# Save session information for reproducibility
save_session <- function(name) {
  init_dirs()
  path <- file.path(PATHS$session, paste0(name, "_session_info.txt"))
  writeLines(capture.output(sessionInfo()), path)
}
