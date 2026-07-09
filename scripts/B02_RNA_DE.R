# B02_RNA_DE.R
# Differential expression analysis with RUVSeq normalization

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
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
})

# Load collapsed counts and metadata
meta_rna <- load_cohort("RNA")

counts_all <- read.table(PATHS$counts_collapsed,
                         header = TRUE, row.names = 1, check.names = FALSE,
                         sep = "\t")
gene_map   <- read.table(PATHS$gene_map, header = TRUE, row.names = 1,
                         sep = "\t")

keep_samples <- intersect(rownames(meta_rna), colnames(counts_all))
counts_rna   <- counts_all[, keep_samples]
meta_rna     <- meta_rna[keep_samples, , drop = FALSE]

stopifnot(identical(colnames(counts_rna), rownames(meta_rna)))
message(sprintf("  %d KS_I | %d control | %d genes",
                sum(meta_rna$status == "KS_I"),
                sum(meta_rna$status == "control"),
                nrow(counts_rna)))

# RUV pass 1: identify empirical control genes
message(sprintf(">>> RUV pass 1 (k = %d)...", K_FACTOR))

dds_p1 <- DESeqDataSetFromMatrix(
  countData = round(as.matrix(counts_rna)),
  colData   = meta_rna,
  design    = ~ batch + sex + status
)
dds_p1 <- dds_p1[rowSums(counts(dds_p1) >= 10) >= 3, ]
dds_p1 <- DESeq(dds_p1, quiet = TRUE)

# Select genes with highest p-values as empirical controls
empirical_genes <- as.data.frame(results(dds_p1)) %>%
  filter(!is.na(pvalue)) %>%
  arrange(desc(pvalue)) %>%
  head(N_EMPIRICAL) %>%
  rownames()

# RUVg normalization using empirical controls
set     <- newSeqExpressionSet(
  as.matrix(counts(dds_p1)),
  phenoData = data.frame(meta_rna, row.names = rownames(meta_rna))
)
set     <- betweenLaneNormalization(set, which = "upper")
set_ruv <- RUVg(set, empirical_genes, k = K_FACTOR)

W_factors <- pData(set_ruv) %>% dplyr::select(starts_with("W_"))
meta_ruv  <- cbind(meta_rna, W_factors)

# RLE plot: before and after RUV normalization
pdf(file.path(PATHS$plots, "B02_RNA_RLE_normalisation.pdf"), width = 10, height = 8)
par(mfrow = c(2, 1))
plotRLE(set,     outline = FALSE, ylim = c(-1, 1),
        col  = as.numeric(meta_rna$status),
        main = "Before RUV")
plotRLE(set_ruv, outline = FALSE, ylim = c(-1, 1),
        col  = as.numeric(meta_rna$status),
        main = paste0("After RUVg (k = ", K_FACTOR, ")"))
dev.off()
message("  -> ", file.path(PATHS$plots, "B02_RNA_RLE_normalisation.pdf"))

# DESeq2 with RUV factors: design = ~ W_1 + W_2 + sex + status
w_terms        <- paste(paste0("W_", seq_len(K_FACTOR)), collapse = " + ")
design_formula <- as.formula(paste0("~ ", w_terms, " + sex + status"))

dds_ruv <- DESeqDataSetFromMatrix(
  countData = counts(dds_p1),
  colData   = meta_ruv,
  design    = design_formula
)
mcols(dds_ruv)$symbol <- gene_map[rownames(dds_ruv), "gene_name"]
dds_ruv <- DESeq(dds_ruv, quiet = TRUE)
saveRDS(dds_ruv, PATHS$dds_rna)

# Shrink log2 fold changes using adaptive shrinkage
res_ruv <- lfcShrink(dds_ruv, coef = "status_KS_I_vs_control", type = "ashr")
res_df  <- as.data.frame(res_ruv) %>%
  rownames_to_column("gene_id") %>%
  mutate(symbol = mcols(dds_ruv)[gene_id, "symbol"]) %>%
  arrange(padj)

write.csv(res_df, PATHS$rna_de, row.names = FALSE)

n_up   <- sum(res_df$padj < P_CUTOFF & res_df$log2FoldChange >  LFC_CUTOFF, na.rm = TRUE)
n_down <- sum(res_df$padj < P_CUTOFF & res_df$log2FoldChange < -LFC_CUTOFF, na.rm = TRUE)

message("  Upregulated   : ", n_up)
message("  Downregulated : ", n_down)
message("  Total DEGs    : ", n_up + n_down)
message("  Threshold     : p.adj < ", P_CUTOFF, " | |Log2FC| > ", LFC_CUTOFF)

# Corrected expression matrix for visualization
vsd     <- vst(dds_ruv, blind = FALSE)
mat_vis <- limma::removeBatchEffect(
  assay(vsd),
  covariates = as.matrix(colData(dds_ruv)[, paste0("W_", seq_len(K_FACTOR)), drop = FALSE]),
  batch2     = dds_ruv$sex,
  design     = model.matrix(~ status, colData(dds_ruv))
)

# PCA on RUV-corrected expression
pca_res    <- prcomp(t(mat_vis))
percentVar <- round(100 * pca_res$sdev^2 / sum(pca_res$sdev^2), 1)
pca_df     <- as.data.frame(pca_res$x) %>%
  rownames_to_column("sample_id") %>%
  left_join(as.data.frame(colData(dds_ruv)) %>% rownames_to_column("sample_id"),
            by = "sample_id")

p_pca <- ggplot(pca_df, aes(PC1, PC2, color = status, shape = sex)) +
  geom_point(size = 4, alpha = 0.9) +
  geom_text_repel(aes(label = sample_id), size = 3, show.legend = FALSE) +
  scale_color_manual(values = COLORS_STATUS) +
  labs(title    = "RNA-seq PCA (RUV corrected)",
       subtitle = sprintf("KS_I: %d | control: %d",
                          sum(pca_df$status == "KS_I"),
                          sum(pca_df$status == "control")),
       x = paste0("PC1: ", percentVar[1], "%"),
       y = paste0("PC2: ", percentVar[2], "%")) +
  theme_bw() + theme(aspect.ratio = 1)

save_plot(p_pca, "B02_RNA_PCA")

# PC3 vs PC4
p_pca34 <- ggplot(pca_df, aes(PC3, PC4, color = status, shape = sex)) +
  geom_point(size = 4, alpha = 0.9) +
  geom_text_repel(aes(label = sample_id), size = 3, show.legend = FALSE) +
  scale_color_manual(values = COLORS_STATUS) +
  labs(title = "RNA-seq PC3 vs PC4 (RUV corrected)",
       x = paste0("PC3: ", percentVar[3], "%"),
       y = paste0("PC4: ", percentVar[4], "%")) +
  theme_bw() + theme(aspect.ratio = 1)

save_plot(p_pca34, "B02_RNA_PCA_PC3_PC4")

# Sample distance heatmap
sampleDist <- as.matrix(dist(t(mat_vis)))
anno_col   <- as.data.frame(colData(dds_ruv)[, c("status", "sex")])

save_heatmap(
  mat         = sampleDist,
  anno_col    = anno_col,
  anno_colors = list(status = COLORS_STATUS, sex = COLORS_SEX),
  title       = "Sample distances (RUV corrected)",
  name        = "B02_RNA_sample_distances",
  show_rownames = TRUE,
  fontsize_row  = 8
)

# Dispersion plot
pdf(file.path(PATHS$plots, "B02_RNA_dispersion.pdf"), width = 7, height = 6)
plotDispEsts(dds_ruv, main = "DESeq2 dispersion estimates")
dev.off()

# MA plot
res_ma <- as.data.frame(results(dds_ruv, name = "status_KS_I_vs_control"))
res_ma$sig <- ifelse(!is.na(res_ma$padj) & res_ma$padj < P_CUTOFF, "FDR<0.05", "ns")

p_ma <- ggplot(res_ma, aes(x = log10(baseMean + 1), y = log2FoldChange,
                            color = sig)) +
  geom_point(alpha = 0.4, size = 0.8) +
  scale_color_manual(values = c("FDR<0.05" = "firebrick", "ns" = "grey70")) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  labs(title = "MA plot (pre-shrinkage)",
       x = "log10(mean normalised count + 1)",
       y = "log2 Fold Change", color = NULL) +
  theme_bw()

save_plot(p_ma, "B02_RNA_MA_plot", width = 7, height = 6)

# Volcano plot with selected gene labels
candidate_genes <- c("VANGL2", "WNT5A", "ROR2", "CELSR1", "PRICKLE1",
                 "ICAM1", "IL6", "ALPL", "PDE10A", "GUCY1A2", "FZD6")
top_genes   <- head(na.omit(res_df$symbol[!is.na(res_df$padj)]), 15)

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
  geom_point(aes(color = Class), alpha = 0.6, size = 1.5) +
  scale_color_manual(values = COLORS_VOLCANO) +
  geom_vline(xintercept = c(-LFC_CUTOFF, LFC_CUTOFF),
             linetype = "dashed", color = "grey50") +
  geom_hline(yintercept = -log10(P_CUTOFF),
             linetype = "dashed", color = "grey50") +
  geom_text_repel(
    data = filter(volcano_df, symbol %in% unique(c(top_genes, candidate_genes))),
    aes(label = symbol), size = 3, max.overlaps = Inf, fontface = "bold"
  ) +
  labs(title    = "RNA-seq: KS_I vs control",
       subtitle = sprintf("%d DEGs (p.adj < %s, |Log2FC| > %s)",
                          n_up + n_down, P_CUTOFF, LFC_CUTOFF),
       x = "Log2 Fold Change", y = "-log10(adjusted p-value)") +
  theme_bw() + theme(aspect.ratio = 1)

save_plot(p_volc, "B02_RNA_volcano")

# Log2 fold change distribution
sig_df <- res_df %>% filter(padj < P_CUTOFF & abs(log2FoldChange) > LFC_CUTOFF)

p_lfc_dist <- ggplot(sig_df, aes(x = log2FoldChange,
                                  fill = log2FoldChange > 0)) +
  geom_histogram(bins = 40, color = "white", linewidth = 0.2) +
  scale_fill_manual(values = c("TRUE" = "firebrick3", "FALSE" = "navy"),
                    labels = c("TRUE" = "Up", "FALSE" = "Down"),
                    name = "Direction") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(title = "LFC distribution of significant DEGs",
       subtitle = sprintf("n = %d (FDR < 0.05, |LFC| > 1)", nrow(sig_df)),
       x = "Log2 Fold Change", y = "Number of genes") +
  theme_bw()

save_plot(p_lfc_dist, "B02_RNA_LFC_distribution", width = 6, height = 5)

# Heatmap of top 100 DEGs
sig_genes <- res_df %>%
  filter(padj < P_CUTOFF & abs(log2FoldChange) > LFC_CUTOFF) %>%
  arrange(padj) %>% head(100)

if (nrow(sig_genes) > 4) {
  mat_heat           <- mat_vis[sig_genes$gene_id, ]
  rownames(mat_heat) <- make.unique(as.character(sig_genes$symbol))

  save_heatmap(
    mat         = mat_heat,
    anno_col    = as.data.frame(colData(dds_ruv)[, c("status", "sex")]),
    anno_colors = list(status = COLORS_STATUS, sex = COLORS_SEX),
    title       = sprintf("Top %d DEGs\n(p.adj < %s, |Log2FC| > %s)",
                          nrow(sig_genes), P_CUTOFF, LFC_CUTOFF),
    name        = "B02_RNA_heatmap_top100"
  )
}

# P-value histogram (QC diagnostic)
p_phist <- ggplot(res_df %>% filter(!is.na(pvalue)),
                  aes(x = pvalue)) +
  geom_histogram(breaks = seq(0, 1, 0.05),
                 fill = "steelblue", color = "white", linewidth = 0.2) +
  labs(title = "P-value histogram (QC)",
       subtitle = "Expected: enrichment near 0 (anti-conservative) for good signal",
       x = "Raw p-value", y = "Number of genes") +
  theme_bw()

save_plot(p_phist, "B02_RNA_pvalue_histogram", width = 6, height = 5)

# RUV factor scatter plot
ruv_df <- as.data.frame(colData(dds_ruv)) %>% rownames_to_column("sample_id")

p_ruv <- ggplot(ruv_df, aes(W_1, W_2, color = status, shape = batch)) +
  geom_point(size = 4) +
  geom_text_repel(aes(label = sample_id), size = 2.8, show.legend = FALSE) +
  scale_color_manual(values = COLORS_STATUS) +
  scale_shape_manual(values = c("84" = 16, "185" = 17)) +
  labs(title = "RUV factors W1 vs W2",
       subtitle = "Should not perfectly align with status — check for confounding",
       x = "W_1", y = "W_2") +
  theme_bw() + theme(aspect.ratio = 1)

save_plot(p_ruv, "B02_RNA_RUV_factors", width = 7, height = 6)

# Background universe for GO enrichment: genes that passed expression
# filtering and were actually tested for DE (not all annotated genes)
universe_ids <- bitr(unique(na.omit(res_df$symbol)), fromType = "SYMBOL",
                     toType = "ENTREZID", OrgDb = org.Hs.eg.db, drop = TRUE)$ENTREZID
saveRDS(universe_ids, file.path(PATHS$rds, "B02_RNA_GO_universe_entrez.rds"))

# GO enrichment (Biological Process) with soft LFC threshold
deg_symbols <- res_df %>%
  filter(padj < P_CUTOFF & abs(log2FoldChange) > LFC_CUTOFF_SOFT) %>%
  pull(symbol) %>% na.omit() %>% unique()

gene_ids <- bitr(deg_symbols, fromType = "SYMBOL", toType = "ENTREZID",
                 OrgDb = org.Hs.eg.db, drop = TRUE)

ego <- enrichGO(
  gene          = gene_ids$ENTREZID,
  universe      = universe_ids,
  OrgDb         = org.Hs.eg.db,
  ont           = GO_ONT,
  pAdjustMethod = "BH",
  pvalueCutoff  = P_CUTOFF,
  readable      = TRUE
)
ego_s <- clusterProfiler::simplify(ego, cutoff = GO_SIMPLIFY,
                                   by = "p.adjust", select_fun = min)

write.csv(as.data.frame(ego_s), PATHS$rna_go, row.names = FALSE)

p_go <- dotplot(ego_s, showCategory = 20) +
  ggtitle("GO Biological Process enrichment",
          subtitle = sprintf("Input: %d DEGs (p.adj < %s, |Log2FC| > %s)",
                             length(deg_symbols), P_CUTOFF, LFC_CUTOFF_SOFT)) +
  theme(plot.title = element_text(size = 11))

save_plot(p_go, "B02_RNA_GO_dotplot", width = 9, height = 10)

# GO enrichment with strict LFC threshold
deg_strict <- res_df %>%
  filter(padj < P_CUTOFF & abs(log2FoldChange) > LFC_CUTOFF) %>%
  pull(symbol) %>% na.omit() %>% unique()

if (length(deg_strict) >= 10) {
  ids_strict <- bitr(deg_strict, fromType = "SYMBOL", toType = "ENTREZID",
                     OrgDb = org.Hs.eg.db, drop = TRUE)
  ego_strict <- enrichGO(gene = ids_strict$ENTREZID, universe = universe_ids,
                         OrgDb = org.Hs.eg.db,
                         ont = GO_ONT, pAdjustMethod = "BH",
                         pvalueCutoff = P_CUTOFF, readable = TRUE)
  if (!is.null(ego_strict) && nrow(ego_strict) > 0) {
    ego_strict_s <- clusterProfiler::simplify(ego_strict, cutoff = GO_SIMPLIFY,
                                              by = "p.adjust", select_fun = min)
    write.csv(as.data.frame(ego_strict_s),
              file.path(PATHS$tables, "B02_RNA_GO_BP_strict.csv"),
              row.names = FALSE)
    p_go_strict <- dotplot(ego_strict_s, showCategory = 20) +
      ggtitle("GO BP enrichment (strict threshold)",
              subtitle = sprintf("Input: %d DEGs (p.adj < %s, |Log2FC| > %s)",
                                 length(deg_strict), P_CUTOFF, LFC_CUTOFF))
    save_plot(p_go_strict, "B02_RNA_GO_dotplot_strict", width = 9, height = 10)
  }
}

save_session("B02")
