# B03_ATAC_DA.R
# Differential accessibility analysis with RUVSeq normalization

source("scripts/00_config.R")
init_dirs()

suppressPackageStartupMessages({
  library(tidyverse)
  library(DESeq2)
  library(RUVSeq)
  library(EDASeq)
  library(limma)
  library(pheatmap)
  library(ggrepel)
  library(DiffBind)
  library(ChIPseeker)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(clusterProfiler)
  library(org.Hs.eg.db)
})

# Load ATAC cohort metadata
meta_atac <- load_cohort("ATAC")

# Extract count matrix from DiffBind object
dbObj    <- readRDS(PATHS$atac_diffbind)
gr_peaks <- dba.peakset(dbObj, bRetrieve = TRUE)

counts_matrix <- as.data.frame(mcols(gr_peaks))
rownames(counts_matrix) <- paste(seqnames(gr_peaks),
                                 start(gr_peaks),
                                 end(gr_peaks), sep = "_")

colnames(counts_matrix) <- dbObj$samples$SampleID

keep_samples <- intersect(rownames(meta_atac), colnames(counts_matrix))
counts_atac  <- counts_matrix[, keep_samples]
meta_atac    <- meta_atac[keep_samples, , drop = FALSE]

stopifnot(identical(colnames(counts_atac), rownames(meta_atac)))
message(sprintf("  Aligned: %d peaks x %d samples", nrow(counts_atac), ncol(counts_atac)))

# RUV pass 1: identify empirical control peaks
message(sprintf(">>> RUV pass 1 (k = %d)...", K_FACTOR))

dds_p1 <- DESeqDataSetFromMatrix(
  countData = round(as.matrix(counts_atac)),
  colData   = meta_atac,
  design    = ~ sex + status
)
dds_p1 <- dds_p1[rowSums(counts(dds_p1) >= 10) >= 3, ]
dds_p1 <- DESeq(dds_p1, quiet = TRUE)

empirical_peaks <- as.data.frame(results(dds_p1)) %>%
  filter(!is.na(pvalue)) %>%
  arrange(desc(pvalue)) %>%
  head(N_EMPIRICAL) %>%
  rownames()

# RUVg normalization
counts_int <- as.matrix(round(counts(dds_p1)))
set        <- newSeqExpressionSet(
  counts_int,
  phenoData = data.frame(meta_atac, row.names = rownames(meta_atac))
)
set        <- betweenLaneNormalization(set, which = "upper")
set_ruv    <- RUVg(set, empirical_peaks, k = K_FACTOR)

W_factors <- pData(set_ruv) %>% dplyr::select(starts_with("W_"))
meta_ruv  <- cbind(meta_atac, W_factors)

# RLE plot
pdf(file.path(PATHS$plots, "B03_ATAC_RLE_normalisation.pdf"), width = 10, height = 8)
par(mfrow = c(2, 1))
plotRLE(set,     outline = FALSE, ylim = c(-1, 1),
        col  = as.numeric(meta_atac$status),
        main = "Before RUV")
plotRLE(set_ruv, outline = FALSE, ylim = c(-1, 1),
        col  = as.numeric(meta_atac$status),
        main = paste0("After RUVg (k = ", K_FACTOR, ")"))
dev.off()
message("  -> ", file.path(PATHS$plots, "B03_ATAC_RLE_normalisation.pdf"))

# DESeq2 with RUV factors: design = ~ W_1 + W_2 + sex + status
w_terms        <- paste(paste0("W_", seq_len(K_FACTOR)), collapse = " + ")
design_formula <- as.formula(paste0("~ ", w_terms, " + sex + status"))

dds_ruv <- DESeqDataSetFromMatrix(
  countData = counts(dds_p1),
  colData   = meta_ruv,
  design    = design_formula
)
dds_ruv <- DESeq(dds_ruv, quiet = TRUE)
saveRDS(dds_ruv, PATHS$dds_atac)

res    <- lfcShrink(dds_ruv, coef = "status_KS_I_vs_control", type = "ashr")
res_df <- as.data.frame(res) %>%
  rownames_to_column("peak_id") %>%
  arrange(padj)

write.csv(res_df, PATHS$atac_da, row.names = FALSE)

n_up   <- sum(res_df$padj < P_CUTOFF & res_df$log2FoldChange >  LFC_CUTOFF, na.rm = TRUE)
n_down <- sum(res_df$padj < P_CUTOFF & res_df$log2FoldChange < -LFC_CUTOFF, na.rm = TRUE)

message("  Opening peaks (UP)  : ", n_up)
message("  Closing peaks (DOWN): ", n_down)
message("  Total DA peaks      : ", n_up + n_down)
message("  UP/DOWN ratio       : ", round(n_up / max(n_down, 1), 2))

# Peak annotation using ChIPseeker
res_split <- res_df %>%
  tidyr::extract(peak_id, into = c("chr", "start", "end"),
                 regex = "^(.+)_([0-9]+)_([0-9]+)$", remove = FALSE) %>%
  mutate(start = as.integer(start), end = as.integer(end))

gr_all <- makeGRangesFromDataFrame(res_split, keep.extra.columns = TRUE)
anno   <- annotatePeak(gr_all,
                       TxDb    = TxDb.Hsapiens.UCSC.hg38.knownGene,
                       annoDb  = "org.Hs.eg.db",
                       verbose = FALSE)
anno_df <- as.data.frame(anno)
write.csv(anno_df, PATHS$atac_annotated, row.names = FALSE)

# Corrected expression matrix for visualization
vsd     <- vst(dds_ruv, blind = FALSE)
mat_vis <- limma::removeBatchEffect(
  assay(vsd),
  covariates = as.matrix(colData(dds_ruv)[, paste0("W_", seq_len(K_FACTOR)), drop = FALSE]),
  batch2     = dds_ruv$sex,
  design     = model.matrix(~ status, colData(dds_ruv))
)

# PCA
pca_res    <- prcomp(t(mat_vis))
percentVar <- round(100 * pca_res$sdev^2 / sum(pca_res$sdev^2), 1)

pca_df <- as.data.frame(pca_res$x) %>%
  rownames_to_column("gdb") %>%
  left_join(as.data.frame(meta_ruv) %>% tibble::rownames_to_column("sample_id"), by = c("gdb" = "sample_id"))

p_pca <- ggplot(pca_df, aes(PC1, PC2, color = status, shape = sex)) +
  geom_point(size = 4, alpha = 0.9) +
  geom_text_repel(aes(label = gdb), size = 3, show.legend = FALSE) +
  scale_color_manual(values = COLORS_STATUS) +
  labs(title    = "ATAC-seq PCA (RUV corrected)",
       subtitle = sprintf("KS_I: %d | control: %d",
                          sum(pca_df$status == "KS_I"),
                          sum(pca_df$status == "control")),
       x = paste0("PC1: ", percentVar[1], "%"),
       y = paste0("PC2: ", percentVar[2], "%")) +
  theme_bw() + theme(aspect.ratio = 1)

save_plot(p_pca, "B03_ATAC_PCA")

# Sample distance heatmap
sampleDist <- as.matrix(dist(t(mat_vis)))
anno_col   <- as.data.frame(colData(dds_ruv)[, c("status", "sex")])

save_heatmap(
  mat         = sampleDist,
  anno_col    = anno_col,
  anno_colors = list(status = COLORS_STATUS, sex = COLORS_SEX),
  title       = "ATAC-seq sample distances (RUV corrected)",
  name        = "B03_ATAC_sample_distances",
  show_rownames = TRUE,
  fontsize_row  = 9
)

# Dispersion plot
pdf(file.path(PATHS$plots, "B03_ATAC_dispersion.pdf"), width = 7, height = 6)
plotDispEsts(dds_ruv, main = "DESeq2 dispersion estimates (ATAC)")
dev.off()
message("  -> ", file.path(PATHS$plots, "B03_ATAC_dispersion.pdf"))

# Volcano plot
volcano_df <- res_df %>%
  mutate(
    Class = case_when(
      padj < P_CUTOFF & log2FoldChange >  LFC_CUTOFF ~ "Upregulated",
      padj < P_CUTOFF & log2FoldChange < -LFC_CUTOFF ~ "Downregulated",
      TRUE ~ "Not Significant"
    ),
    neg_log10_padj = pmin(-log10(padj), 300)
  )

p_volc <- ggplot(volcano_df, aes(log2FoldChange, neg_log10_padj)) +
  geom_point(aes(color = Class), alpha = 0.6, size = 1.0) +
  scale_color_manual(values = COLORS_VOLCANO) +
  geom_vline(xintercept = c(-LFC_CUTOFF, LFC_CUTOFF),
             linetype = "dashed", color = "grey50") +
  geom_hline(yintercept = -log10(P_CUTOFF),
             linetype = "dashed", color = "grey50") +
  labs(title    = "ATAC-seq: KS_I vs control",
       subtitle = sprintf("%d DA peaks (p.adj < %s, |Log2FC| > %s)",
                          n_up + n_down, P_CUTOFF, LFC_CUTOFF),
       x = "Log2 Fold Change", y = "-log10(adjusted p-value)") +
  theme_bw() + theme(aspect.ratio = 1)

save_plot(p_volc, "B03_ATAC_volcano")

# MA plot
res_ma <- as.data.frame(results(dds_ruv, name = "status_KS_I_vs_control"))
res_ma$sig <- ifelse(!is.na(res_ma$padj) & res_ma$padj < P_CUTOFF, "FDR<0.05", "ns")

p_ma <- ggplot(res_ma, aes(x = log10(baseMean + 1), y = log2FoldChange,
                            color = sig)) +
  geom_point(alpha = 0.3, size = 0.6) +
  scale_color_manual(values = c("FDR<0.05" = "firebrick", "ns" = "grey70")) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  labs(title = "ATAC-seq MA plot (pre-shrinkage)",
       x = "log10(mean normalised count + 1)",
       y = "log2 Fold Change", color = NULL) +
  theme_bw()

save_plot(p_ma, "B03_ATAC_MA_plot", width = 7, height = 6)

# Accessibility asymmetry barplot
df_asym <- data.frame(
  Direction = c("Opening (UP)", "Closing (DOWN)"),
  Count     = c(n_up, n_down)
)

p_asym <- ggplot(df_asym, aes(Direction, Count, fill = Direction)) +
  geom_bar(stat = "identity", width = 0.6, color = "black") +
  scale_fill_manual(values = c("Opening (UP)"  = "firebrick3",
                               "Closing (DOWN)" = "navy")) +
  geom_text(aes(label = Count), vjust = -0.5, size = 5, fontface = "bold") +
  labs(title    = "Chromatin accessibility asymmetry",
       subtitle = sprintf("UP/DOWN = %.2fx", n_up / max(n_down, 1)),
       y = "Number of DA peaks", x = "") +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 12, face = "bold"))

save_plot(p_asym, "B03_ATAC_asymmetry", width = 5, height = 6)

# Peak annotation comparison
simplify_feature <- function(x) {
  case_when(
    grepl("Promoter",   x) ~ "Promoter",
    grepl("Exon",       x) ~ "Exon",
    grepl("Intron",     x) ~ "Intron",
    grepl("Downstream", x) ~ "Downstream",
    grepl("Distal",     x) ~ "Distal Intergenic",
    TRUE                   ~ "Other"
  )
}

feature_colors <- c(
  "Promoter"          = "#e74c3c",
  "Intron"            = "#3498db",
  "Distal Intergenic" = "#27ae60",
  "Exon"              = "#f39c12",
  "Downstream"        = "#8e44ad",
  "Other"             = "grey70"
)

sig_peak_ids <- res_df %>%
  filter(padj < P_CUTOFF & abs(log2FoldChange) > LFC_CUTOFF) %>%
  pull(peak_id)

anno_summary <- bind_rows(
  anno_df %>%
    mutate(feature = simplify_feature(annotation), set = "All tested peaks") %>%
    dplyr::count(set, feature),
  anno_df %>%
    filter(peak_id %in% sig_peak_ids) %>%
    mutate(feature = simplify_feature(annotation), set = "DA peaks (FDR<0.05, |LFC|>1)") %>%
    dplyr::count(set, feature)
) %>%
  group_by(set) %>%
  mutate(pct     = 100 * n / sum(n),
         feature = factor(feature, levels = names(feature_colors))) %>%
  ungroup()

p_anno <- ggplot(anno_summary, aes(x = set, y = pct, fill = feature)) +
  geom_bar(stat = "identity", width = 0.5, colour = "white", linewidth = 0.4) +
  geom_text(aes(label = ifelse(pct > 5, paste0(round(pct, 1), "%"), "")),
            position = position_stack(vjust = 0.5),
            size = 3, colour = "white", fontface = "bold") +
  scale_fill_manual(values = feature_colors, name = "Genomic feature") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 101)) +
  labs(title    = "Genomic annotation of ATAC-seq peaks",
       subtitle = sprintf("Tested: %s | DA: %d",
                          format(nrow(res_df), big.mark = ","),
                          length(sig_peak_ids)),
       x = NULL, y = "Peaks (%)") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 9.5))

save_plot(p_anno, "B03_ATAC_peak_annotation", width = 5, height = 5)

# P-value histogram
p_phist <- ggplot(res_df %>% filter(!is.na(pvalue)), aes(x = pvalue)) +
  geom_histogram(breaks = seq(0, 1, 0.05),
                 fill = "steelblue", color = "white", linewidth = 0.2) +
  labs(title = "ATAC p-value histogram (QC)",
       x = "Raw p-value", y = "Number of peaks") +
  theme_bw()

save_plot(p_phist, "B03_ATAC_pvalue_histogram", width = 6, height = 5)

# RUV factors
ruv_df <- as.data.frame(colData(dds_ruv)) %>% rownames_to_column("sample_id")

p_ruv <- ggplot(ruv_df, aes(W_1, W_2, color = status, shape = sex)) +
  geom_point(size = 4) +
  geom_text_repel(aes(label = sample_id), size = 2.8, show.legend = FALSE) +
  scale_color_manual(values = COLORS_STATUS) +
  labs(title = "RUV factors W1 vs W2 (ATAC)",
       x = "W_1", y = "W_2") +
  theme_bw() + theme(aspect.ratio = 1)

save_plot(p_ruv, "B03_ATAC_RUV_factors", width = 7, height = 6)

# GO enrichment for opening vs closing peaks
run_atac_go <- function(peaks_sub, label) {
  if (nrow(peaks_sub) < 20) {
    message("  Skipping GO [", label, "] — n=", nrow(peaks_sub), " peaks")
    return(NULL)
  }
  genes <- anno_df %>%
    filter(peak_id %in% peaks_sub$peak_id) %>%
    pull(geneId) %>% unique() %>% na.omit()
  ego <- enrichGO(gene = genes, OrgDb = org.Hs.eg.db, ont = GO_ONT,
                  pAdjustMethod = "BH", readable = TRUE)
  if (is.null(ego) || nrow(ego) == 0) return(NULL)
  ego_s <- clusterProfiler::simplify(ego, cutoff = GO_SIMPLIFY,
                                     by = "p.adjust", select_fun = min)
  write.csv(as.data.frame(ego_s),
            file.path(PATHS$tables, paste0("B03_ATAC_GO_", label, ".csv")),
            row.names = FALSE)
  p <- dotplot(ego_s, showCategory = 15) +
    ggtitle(paste0("GO BP: ", label, " peaks (n=", length(genes), " genes)"))
  save_plot(p, paste0("B03_ATAC_GO_", label), width = 9, height = 9)
  return(as.data.frame(ego_s))
}

peaks_up   <- res_df %>% filter(padj < P_CUTOFF & log2FoldChange >  LFC_CUTOFF)
peaks_down <- res_df %>% filter(padj < P_CUTOFF & log2FoldChange < -LFC_CUTOFF)

run_atac_go(peaks_up,   "opening")
run_atac_go(peaks_down, "closing")

# LFC threshold sensitivity analysis
for (lfc in c(1.0, 0.75, 0.585, 0.5)) {
  nu <- sum(res_df$padj < 0.05 & res_df$log2FoldChange >  lfc, na.rm = TRUE)
  nd <- sum(res_df$padj < 0.05 & res_df$log2FoldChange < -lfc, na.rm = TRUE)
  message(sprintf("  LFC > %.3f : UP=%d | DOWN=%d | Total=%d", lfc, nu, nd, nu + nd))
}

save_session("B03")
