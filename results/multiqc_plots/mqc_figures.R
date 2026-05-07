#!/usr/bin/env Rscript
# mqc_figures.R extract MultiQC data and plot QC figures for Scientific Data
# Dependencies: lzstring, jsonlite, ggplot2, dplyr, tidyr, patchwork
# Install once: install.packages(c("lzstring","jsonlite","ggplot2","dplyr","tidyr","patchwork"))

library(lzstring)
library(jsonlite)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

# ── Config ────────────────────────────────────────────────────────────────────
ATAC_HTML <- "ATACmultiqc_report.html"
RNA_HTML  <- "RNAfibro multiqc_report.html"
dir.create("figures/ATAC", recursive = TRUE, showWarnings = FALSE)
dir.create("figures/RNA",  recursive = TRUE, showWarnings = FALSE)

CTRL_COL <- "#4393C3"
KS_COL   <- "#D6604D"
CB_CATS  <- c("#437bb1","#b1084c","#4CAF50","#E87722","#7b3294","#888888","#f4a582","#92c5de")

THEME <- theme_classic() +
  theme(axis.text.y = element_text(size = 7),
        plot.title  = element_text(size = 9, face = "bold"),
        legend.text = element_text(size = 7),
        legend.key.size = unit(0.4, "cm"))

save_fig <- function(p, path, name, w = 10, h = 4) {
  ggsave(file.path(path, paste0(name, ".png")), p, width = w, height = h, dpi = 300)
  ggsave(file.path(path, paste0(name, ".svg")), p, width = w, height = h)
  message("  saved  ", file.path(path, name))
}

# ── Load MultiQC ──────────────────────────────────────────────────────────────
load_mqc <- function(html_path) {
  html    <- readChar(html_path, file.info(html_path)$size, useBytes = FALSE)
  scripts <- regmatches(html, gregexpr("(?s)<script[^>]*>(.*?)</script>", html, perl = TRUE))[[1]]
  raw     <- gsub("(?s)^<script[^>]*>|</script>$", "", scripts[1], perl = TRUE)
  fromJSON(lzstring::decompressFromBase64(trimws(raw)), simplifyVector = FALSE)
}

# ── Helpers ───────────────────────────────────────────────────────────────────
# bar_graph: plot$samples[[1]] = sample names; plot$datasets[[1]] = list of {name, data}
bar_to_df <- function(plot, idx = 1) {
  samples <- unlist(plot$samples[[idx]])
  cats    <- plot$datasets[[idx]]
  df <- as.data.frame(sapply(cats, function(x) unlist(x[["data"]])))
  colnames(df) <- sapply(cats, `[[`, "name")
  df$sample <- samples
  df
}

# xy_line: plot$datasets[[1]] = list of {name, data:[[x,y],...]}
xyline_to_df <- function(plot, idx = 1) {
  lapply(plot$datasets[[idx]], function(s) {
    pts <- s[["data"]]
    m   <- do.call(rbind, lapply(pts, function(p) as.numeric(unlist(p))))
    data.frame(sample = s[["name"]], x = m[, 1], y = m[, 2])
  }) |> bind_rows()
}

# scatter: plot$datasets[[1]] = list of {x, y, name}
scatter_to_df <- function(plot, idx = 1) {
  pts <- plot$datasets[[idx]]
  data.frame(x    = sapply(pts, function(p) as.numeric(p[["x"]])),
             y    = sapply(pts, function(p) as.numeric(p[["y"]])),
             name = sapply(pts, `[[`, "name"))
}

# sparse diagonal format (FRiP, peak count)
sparse_to_df <- function(plot, idx = 1) {
  samples <- unlist(plot$samples[[idx]])
  cats    <- plot$datasets[[idx]]
  vals <- sapply(seq_along(samples), function(i) {
    cat <- Filter(function(x) x[["name"]] == samples[i], cats)
    if (length(cat)) as.numeric(cat[[1]][["data"]][[i]]) else NA
  })
  data.frame(sample = samples, value = vals)
}

# Lane merging: sum across lanes/read-pairs then groupby base name
strip_lane_rp <- function(s) sub("_[12]$", "", sub("_T\\d+$", "", s))  # ATAC FastQC
strip_lane    <- function(s) sub("_T\\d+$", "", s)                      # ATAC Samtools
strip_rp      <- function(s) sub("_[12]$", "", s)                       # RNA FastQC

merge_bar <- function(df, key_fn) {
  df$sample <- key_fn(df$sample)
  df |> group_by(sample) |> summarise(across(where(is.numeric), sum), .groups = "drop")
}

merge_xy <- function(df, key_fn) {
  df$sample <- key_fn(df$sample)
  df |> group_by(sample, x) |> summarise(y = mean(y), .groups = "drop")
}

cond_col <- function(s) ifelse(grepl("ctrl", s), CTRL_COL, KS_COL)

# ── Plot functions ─────────────────────────────────────────────────────────────
stacked_barh <- function(df, title, xlab, pct = FALSE) {
  long <- df |> pivot_longer(-sample, names_to = "cat", values_to = "val")
  if (pct) long <- long |> group_by(sample) |> mutate(val = val / sum(val) * 100) |> ungroup()
  ggplot(long, aes(y = sample, x = val, fill = cat)) +
    geom_col() +
    scale_fill_manual(values = setNames(CB_CATS[seq_along(unique(long$cat))], unique(long$cat))) +
    labs(title = title, x = if (pct) "%" else xlab, y = NULL, fill = NULL) +
    THEME
}

line_p <- function(df, title, xlab, ylab) {
  ggplot(df, aes(x = x, y = y, colour = sample)) +
    geom_line(linewidth = 0.7) +
    geom_hline(yintercept = 28, linetype = "dashed", colour = "orange", linewidth = 0.6) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    coord_cartesian(ylim = c(0, NA)) +
    labs(title = title, x = xlab, y = ylab, colour = NULL) +
    guides(colour = guide_legend(ncol = 2, keyheight = 0.5)) +
    THEME + theme(legend.text = element_text(size = 5.5))
}

pca_p <- function(df, title, xlab = "PC1", ylab = "PC2") {
  df$cond <- cond_col(df$name)
  ggplot(df, aes(x = x, y = y, colour = cond, label = name)) +
    geom_point(size = 2.5) +
    geom_text(size = 2, hjust = -0.1, vjust = 0.5) +
    scale_colour_identity(guide = "legend",
                          labels = c(CTRL_COL = "ctrl", KS_COL = "KS"),
                          breaks = c(CTRL_COL, KS_COL)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey70", linewidth = 0.4) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey70", linewidth = 0.4) +
    labs(title = title, x = xlab, y = ylab, colour = NULL) +
    THEME
}

single_barh <- function(df, title, xlab, vline = NULL) {
  df$cond <- cond_col(df$sample)
  p <- ggplot(df, aes(y = sample, x = value, fill = cond)) +
    geom_col() +
    scale_fill_identity(guide = "legend",
                        labels = c(CTRL_COL = "ctrl", KS_COL = "KS"),
                        breaks = c(CTRL_COL, KS_COL)) +
    labs(title = title, x = xlab, y = NULL, fill = NULL) +
    THEME
  if (!is.null(vline))
    p <- p + geom_vline(xintercept = vline, linetype = "dashed", colour = "red", linewidth = 0.8)
  p
}

# ══════════════════════════════════════════════════════════════════════════════
# ATAC-seq
# ══════════════════════════════════════════════════════════════════════════════
message("\n── ATAC-seq ──")
atac <- load_mqc(ATAC_HTML)

# FastQC per-base quality post-trim (lanes + read-pairs merged)
df <- xyline_to_df(atac[["fastqc_per_base_sequence_quality_plot-2-1"]]) |>
  merge_xy(strip_lane_rp)
save_fig(line_p(df, "FastQC – Per-base quality (post-trim)", "Position (bp)", "Phred score"),
         "figures/ATAC", "A1_FastQC_quality_posttrim")

# Cutadapt filtered reads (lanes merged)
df <- bar_to_df(atac[["cutadapt_filtered_reads_plot-1"]]) |> merge_bar(strip_lane)
save_fig(stacked_barh(df, "Cutadapt – Filtered reads", "Reads"),
         "figures/ATAC", "A2_Cutadapt")

# Samtools alignment post-dedup
df <- bar_to_df(atac[["samtools_alignment_plot-3-1"]])
save_fig(stacked_barh(df, "Samtools – Alignment (post-dedup)", "Reads"),
         "figures/ATAC", "A3_Samtools_dedup")

# Picard insert size — normalised to % and smoothed (loess)
df_ins <- xyline_to_df(atac[["picard_insert_size-1"]]) |>
  filter(x >= 0, x <= 800) |>
  group_by(sample) |>
  mutate(y = y / sum(y) * 100) |>
  ungroup()
save_fig(
  ggplot(df_ins, aes(x = x, y = y, colour = sample)) +
    geom_smooth(method = "loess", span = 0.08, se = FALSE, linewidth = 0.8) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    coord_cartesian(ylim = c(0, NA)) +
    labs(title = "Picard – Insert size distribution (%)",
         x = "Insert size (bp)", y = "% of reads", colour = NULL) +
    guides(colour = guide_legend(ncol = 2, keyheight = 0.5)) +
    THEME + theme(legend.text = element_text(size = 5.5)),
  "figures/ATAC", "A4_Picard_insert_size"
)

# DeepTools fingerprint
df <- xyline_to_df(atac[["deeptools_fingerprint_plot-1"]])
save_fig(line_p(df, "DeepTools – Fingerprint", "Fraction of genome", "Fraction of reads"),
         "figures/ATAC", "A5_Fingerprint")

# FRiP score
df <- sparse_to_df(atac[["mlib_frip_score-plot-1"]])
save_fig(single_barh(df, "FRiP score", "FRiP", vline = 0.2),
         "figures/ATAC", "A6_FRiP", w = 7)

# Peak annotation
df <- bar_to_df(atac[["mlib_peak_annotation-plot-1"]])
save_fig(stacked_barh(df, "Peak annotation (%)", "%", pct = TRUE),
         "figures/ATAC", "A7_Peak_annotation")

# DESeq2 PCA
df  <- scatter_to_df(atac[["mlib_deseq2_pca_1-plot"]])
cfg <- atac[["mlib_deseq2_pca_1-plot"]]$config
save_fig(pca_p(df, "DESeq2 PCA – ATAC-seq",
               cfg$xlab %||% "PC1", cfg$ylab %||% "PC2"),
         "figures/ATAC", "A8_DESeq2_PCA", w = 7, h = 6)

# ══════════════════════════════════════════════════════════════════════════════
# RNA-seq
# ══════════════════════════════════════════════════════════════════════════════
message("\n── RNA-seq ──")
rna <- load_mqc(RNA_HTML)

# FastQC per-base quality post-trim
df <- xyline_to_df(rna[["fastqc_per_base_sequence_quality_plot-2-1"]]) |>
  merge_xy(strip_rp)
save_fig(line_p(df, "FastQC – Per-base quality (post-trim)", "Position (bp)", "Phred score"),
         "figures/RNA", "R1_FastQC_quality_posttrim")

# Cutadapt
df <- bar_to_df(rna[["cutadapt_filtered_reads_plot-1"]])
save_fig(stacked_barh(df, "Cutadapt – Filtered reads", "Reads"),
         "figures/RNA", "R2_Cutadapt")

# STAR alignment
df <- bar_to_df(rna[["star_alignment_plot-1"]])
save_fig(stacked_barh(df, "STAR – Alignment summary", "Reads"),
         "figures/RNA", "R3_STAR_alignment")

# RSeQC read distribution
df <- bar_to_df(rna[["rseqc_read_distribution_plot-1"]])
save_fig(stacked_barh(df, "RSeQC – Read distribution (%)", "%", pct = TRUE),
         "figures/RNA", "R4_RSeQC_distribution")

# RSeQC infer experiment
df <- bar_to_df(rna[["rseqc_infer_experiment_plot-1"]])
save_fig(stacked_barh(df, "RSeQC – Strandedness (%)", "%", pct = TRUE),
         "figures/RNA", "R5_RSeQC_strandedness")

# DESeq2 PCA
df  <- scatter_to_df(rna[["star_salmon_deseq2_pca-plot"]])
cfg <- rna[["star_salmon_deseq2_pca-plot"]]$config
save_fig(pca_p(df, "DESeq2 PCA – RNA-seq",
               cfg$xlab %||% "PC1", cfg$ylab %||% "PC2"),
         "figures/RNA", "R6_DESeq2_PCA", w = 7, h = 6)

message("\nDone. 14 figures in figures/ATAC/ and figures/RNA/")