#!/usr/bin/env Rscript
# DiffBind differential accessibility analysis
# GRCh38 | broad peaks | two normalization strategies

suppressPackageStartupMessages({
  library(DiffBind)
  library(rtracklayer)
  library(BiocParallel)
  library(ggplot2)
  library(ChIPseeker)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(org.Hs.eg.db)
  library(GenomicRanges)
  library(ggrepel)
  library(pheatmap)
  library(RColorBrewer)
})

register(MulticoreParam(workers = 16))

# Project directory structure
PROJECT   <- file.path(Sys.getenv("HOME"), "analysis/KABUKI_FIBRO_PROJECT")
ATAC_RES  <- file.path(PROJECT, "pipeline/ATAC/results")
BAM_DIR   <- file.path(ATAC_RES, "bwa/merged_library")
PEAK_DIR  <- file.path(BAM_DIR,  "macs2/broad_peak")
OUT_DIR   <- file.path(PROJECT,  "results/ATAC/DiffBind")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT_DIR, "plots"),  showWarnings = FALSE)
dir.create(file.path(OUT_DIR, "tables"), showWarnings = FALSE)
dir.create(file.path(OUT_DIR, "beds"),   showWarnings = FALSE)

# Sample sheet
meta <- data.frame(
  SampleID  = c("GDB339_ctrl", "GDB418_ctrl", "GDB809_ctrl", "GDB813_ctrl",
                "GDB374_KS",   "GDB380_KS",   "GDB425_KS",
                "GDB447_KS",   "GDB587_KS",   "GDB598_KS"),
  Condition = c(rep("control", 4), rep("KS", 6)),
  Replicate = 1,
  stringsAsFactors = FALSE
)

meta$bamReads   <- file.path(BAM_DIR,  paste0(meta$SampleID, "_REP1.mLb.clN.sorted.bam"))
meta$Peaks      <- file.path(PEAK_DIR, paste0(meta$SampleID, "_REP1.mLb.clN_peaks.broadPeak"))
meta$PeakCaller <- "bed"
meta$Tissue     <- "fibroblast"
meta$Factor     <- "ATAC"

write.csv(meta, file.path(OUT_DIR, "tables/diffbind_samplesheet.csv"),
          row.names = FALSE, quote = FALSE)

# Verify input files exist
missing_bam  <- meta$bamReads[!file.exists(meta$bamReads)]
missing_peak <- meta$Peaks[!file.exists(meta$Peaks)]
if (length(missing_bam) > 0 | length(missing_peak) > 0) {
  stop("Missing input files. Check paths.")
}

# Load or count peaks (checkpoint for resuming)
counted_rds <- file.path(OUT_DIR, "diffbind_counted.rds")

if (file.exists(counted_rds)) {
  atac <- readRDS(counted_rds)
} else {
  atac <- dba(sampleSheet = meta)
  atac$config$cores <- 16
  atac <- dba.count(atac, bUseSummarizeOverlaps = TRUE)
  saveRDS(atac, counted_rds)
}
atac$config$cores <- 16

# QC metrics
info <- as.data.frame(dba.show(atac))

qc <- data.frame(
  SampleID  = info$ID,
  Condition = info$Condition,
  Reads     = info$Reads,
  FRiP      = round(info$FRiP, 4),
  Peaks     = info$Intervals,
  stringsAsFactors = FALSE
)

write.csv(qc, file.path(OUT_DIR, "tables/QC_metrics.csv"),
          row.names = FALSE, quote = FALSE)

# FRiP barplot
p_frip <- ggplot(qc, aes(x = reorder(SampleID, FRiP), y = FRiP, fill = Condition)) +
  geom_col() +
  geom_hline(yintercept = 0.2, linetype = "dashed", color = "grey40") +
  annotate("text", x = 1, y = 0.21, label = "FRiP = 0.2", hjust = 0,
           size = 3, color = "grey40") +
  scale_fill_manual(values = c("control" = "#4393C3", "KS" = "#D6604D")) +
  coord_flip() +
  labs(title = "Fraction of Reads in Peaks",
       x = NULL, y = "FRiP", fill = "Condition") +
  theme_bw(base_size = 12)
ggsave(file.path(OUT_DIR, "plots/QC_FRiP.pdf"), p_frip, width = 7, height = 5)

# Library size barplot
p_lib <- ggplot(qc, aes(x = reorder(SampleID, Reads), y = Reads / 1e6, fill = Condition)) +
  geom_col() +
  scale_fill_manual(values = c("control" = "#4393C3", "KS" = "#D6604D")) +
  coord_flip() +
  labs(title = "Mapped Read Counts",
       x = NULL, y = "Reads (millions)", fill = "Condition") +
  theme_bw(base_size = 12)
ggsave(file.path(OUT_DIR, "plots/QC_LibrarySize.pdf"), p_lib, width = 7, height = 5)

# Sample correlation heatmap
pdf(file.path(OUT_DIR, "plots/QC_CorrelationHeatmap.pdf"), width = 8, height = 7)
dba.plotHeatmap(atac, ColAttributes = DBA_CONDITION,
                colScheme = "RdYlBu", main = "Sample correlation (all peaks)")
dev.off()

# PCA on all peaks
pdf(file.path(OUT_DIR, "plots/QC_PCA_all_peaks.pdf"), width = 7, height = 6)
dba.plotPCA(atac, attributes = DBA_CONDITION, label = DBA_ID)
dev.off()

# Differential analysis function (runs one normalization strategy)
run_branch <- function(dba_obj, norm_label, out_subdir) {

  dir.create(out_subdir, recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(out_subdir, "plots"),  showWarnings = FALSE)
  dir.create(file.path(out_subdir, "tables"), showWarnings = FALSE)
  dir.create(file.path(out_subdir, "beds"),   showWarnings = FALSE)

  # Apply normalization
  if (norm_label == "Background") {
    dba_obj <- dba.normalize(dba_obj, background = TRUE)
  } else {
    dba_obj <- dba.normalize(dba_obj)
  }

  # Save normalization factors
  norm_info <- dba.normalize(dba_obj, bRetrieve = TRUE)
  norm_df <- data.frame(
    SampleID   = meta$SampleID,
    NormFactor = norm_info$norm.factors
  )
  write.csv(norm_df,
            file.path(out_subdir, "tables/normalization_factors.csv"),
            row.names = FALSE, quote = FALSE)

  # Set contrast and run DESeq2
  dba_obj <- dba.contrast(dba_obj, contrast = c("Condition", "KS", "control"))
  dba_obj <- dba.analyze(dba_obj, method = DBA_DESEQ2)

  saveRDS(dba_obj, file.path(out_subdir, "diffbind_analyzed.rds"))

  # Export results
  res_all <- dba.report(dba_obj, contrast = 1, th = 1,
                         bUsePval = FALSE, bCounts = TRUE)
  res_sig <- dba.report(dba_obj, contrast = 1, th = 0.05,
                         bUsePval = FALSE, bCounts = TRUE)

  df_all <- as.data.frame(res_all)
  df_sig <- as.data.frame(res_sig)

  write.csv(df_all, file.path(out_subdir, "tables/DA_all_regions.csv"),
            row.names = FALSE, quote = FALSE)
  write.csv(df_sig, file.path(out_subdir, "tables/DA_significant_FDR05.csv"),
            row.names = FALSE, quote = FALSE)

  # Export BED files for gained/lost peaks
  if (nrow(df_sig) > 0) {
    gained <- res_sig[res_sig$Fold > 0]
    lost   <- res_sig[res_sig$Fold < 0]
    if (length(gained) > 0) export.bed(gained, file.path(out_subdir, "beds/KS_GAINED.bed"))
    if (length(lost)   > 0) export.bed(lost,   file.path(out_subdir, "beds/KS_LOST.bed"))
  }

  # Volcano plot
  df_all$sig <- "ns"
  df_all$sig[df_all$FDR < 0.05 & df_all$Fold > 0] <- "Gained"
  df_all$sig[df_all$FDR < 0.05 & df_all$Fold < 0] <- "Lost"
  df_all$sig <- factor(df_all$sig, levels = c("Gained", "Lost", "ns"))
  df_all$neglog10fdr <- pmin(-log10(df_all$FDR), 50)

  n_gained <- sum(df_all$sig == "Gained")
  n_lost   <- sum(df_all$sig == "Lost")

  p_vol <- ggplot(df_all, aes(x = Fold, y = neglog10fdr, color = sig)) +
    geom_point(alpha = 0.5, size = 1.0) +
    scale_color_manual(values = c("Gained" = "#D6604D", "Lost" = "#4393C3", "ns" = "grey75"),
                       labels = c(paste0("Gained (n=", n_gained, ")"),
                                  paste0("Lost (n=",   n_lost,   ")"), "ns")) +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", linewidth = 0.4) +
    labs(title = paste0("KS vs Control (", norm_label, " normalization)"),
         x = "Fold change (log2)", y = "-log10(FDR)", color = NULL) +
    theme_bw(base_size = 12) +
    theme(legend.position = "top")
  ggsave(file.path(out_subdir, "plots/Volcano.pdf"), p_vol, width = 7, height = 6)

  # MA plot
  conc_cols <- grep("^Conc_", colnames(df_all), value = TRUE)
  df_all$meanConc <- if (length(conc_cols) >= 2) rowMeans(df_all[, conc_cols]) else df_all$Conc

  p_ma <- ggplot(df_all, aes(x = meanConc, y = Fold, color = sig)) +
    geom_point(alpha = 0.5, size = 1.0) +
    scale_color_manual(values = c("Gained" = "#D6604D", "Lost" = "#4393C3", "ns" = "grey75")) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4) +
    labs(title = paste0("MA plot (", norm_label, " normalization)"),
         x = "Mean concentration (log2)", y = "Fold change (log2)", color = NULL) +
    theme_bw(base_size = 12) +
    theme(legend.position = "top")
  ggsave(file.path(out_subdir, "plots/MA.pdf"), p_ma, width = 7, height = 6)

  # PCA on DA regions
  pdf(file.path(out_subdir, "plots/PCA_DA_regions.pdf"), width = 7, height = 6)
  dba.plotPCA(dba_obj, contrast = 1, attributes = DBA_CONDITION,
              label = DBA_ID)
  dev.off()

  # Correlation heatmap on DA regions
  pdf(file.path(out_subdir, "plots/CorrelationHeatmap_DA.pdf"), width = 8, height = 7)
  dba.plotHeatmap(dba_obj, contrast = 1, ColAttributes = DBA_CONDITION,
                  colScheme = "RdYlBu",
                  main = paste0("Correlation (DA regions, ", norm_label, ")"))
  dev.off()

  # Genomic annotation using ChIPseeker
  if (nrow(df_sig) > 0) {
    txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
    beds_to_annotate <- list()
    gained_bed <- file.path(out_subdir, "beds/KS_GAINED.bed")
    lost_bed   <- file.path(out_subdir, "beds/KS_LOST.bed")
    if (file.exists(gained_bed)) beds_to_annotate[["KS Gained"]] <- import.bed(gained_bed)
    if (file.exists(lost_bed))   beds_to_annotate[["KS Lost"]]   <- import.bed(lost_bed)

    anno_list <- lapply(beds_to_annotate, function(gr) {
      annotatePeak(gr, tssRegion = c(-3000, 3000),
                   TxDb = txdb, annoDb = "org.Hs.eg.db", verbose = FALSE)
    })

    # Annotation barplot
    pdf(file.path(out_subdir, "plots/Annotation_BarChart.pdf"), width = 8, height = 5)
    plotAnnoBar(anno_list, title = paste0("Genomic annotation (", norm_label, ")"))
    dev.off()

    # Distance to TSS distribution
    pdf(file.path(out_subdir, "plots/Annotation_DistToTSS.pdf"), width = 8, height = 5)
    plotDistToTSS(anno_list, title = paste0("Distance to TSS (", norm_label, ")"))
    dev.off()

    # Export annotation tables
    for (nm in names(anno_list)) {
      write.csv(as.data.frame(anno_list[[nm]]),
                file.path(out_subdir, paste0("tables/Annotation_", gsub(" ", "_", nm), ".csv")),
                row.names = FALSE, quote = FALSE)
    }
  }

  invisible(dba_obj)
}

# Run both normalization strategies
run_branch(atac, "Default",    file.path(OUT_DIR, "Default_norm"))
run_branch(atac, "Background", file.path(OUT_DIR, "Background_norm"))

# Save session info
sink(file.path(OUT_DIR, "session_info.txt"))
sessionInfo()
sink()
