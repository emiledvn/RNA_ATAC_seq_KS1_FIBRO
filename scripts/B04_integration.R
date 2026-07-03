# B04_integration.R
# Integration of RNA-seq and ATAC-seq results

source("scripts/00_config.R")
init_dirs()

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggrepel)
  library(pheatmap)
  library(DESeq2)
  library(limma)
  library(grid)
  library(gridExtra)
  library(clusterProfiler)
  library(org.Hs.eg.db)
})

# ── 1. Load data ──────────────────────────────────────────────────────────────

res_rna <- read.csv(PATHS$rna_de) %>%
  dplyr::select(gene_id, symbol, log2FoldChange, padj) %>%
  dplyr::rename(logFC_RNA = log2FoldChange, padj_RNA = padj) %>%
  filter(!is.na(gene_id)) %>%
  mutate(gene_id = sub("\\..*$", "", gene_id)) %>%
  group_by(gene_id) %>%
  slice_min(padj_RNA, n = 1, with_ties = FALSE) %>%
  ungroup()

# All peaks with an ENSEMBL annotation; symbol may be NA
res_atac <- read.csv(PATHS$atac_annotated) %>%
  dplyr::select(peak_id, ENSEMBL, SYMBOL, log2FoldChange, padj, annotation) %>%
  dplyr::rename(
    ensembl_atac = ENSEMBL,
    symbol_atac  = SYMBOL,
    logFC_ATAC   = log2FoldChange,
    padj_ATAC    = padj
  ) %>%
  filter(!is.na(ensembl_atac))

# ── 2. Merge and classify ─────────────────────────────────────────────────────

merged_df <- inner_join(
  res_rna,
  res_atac,
  by = c("gene_id" = "ensembl_atac")
) %>%
  mutate(
    label = case_when(
      !is.na(symbol) & symbol != "" ~ symbol,
      !is.na(symbol_atac) & symbol_atac != "" ~ symbol_atac,
      TRUE ~ gene_id
    ),
    sig_RNA  = padj_RNA  < P_CUTOFF & abs(logFC_RNA)  > LFC_CUTOFF_SOFT,
    sig_ATAC = padj_ATAC < P_CUTOFF & abs(logFC_ATAC) > LFC_CUTOFF_SOFT,
    Direction = case_when(
      sig_RNA & sig_ATAC & logFC_RNA > 0 & logFC_ATAC > 0 ~ "UP-Opening",
      sig_RNA & sig_ATAC & logFC_RNA < 0 & logFC_ATAC < 0 ~ "DOWN-Closing",
      sig_RNA & sig_ATAC & logFC_RNA > 0 & logFC_ATAC < 0 ~ "UP-Closing",
      sig_RNA & sig_ATAC & logFC_RNA < 0 & logFC_ATAC > 0 ~ "DOWN-Opening",
      TRUE ~ "Not significant"
    )
  )

# ── 3. DEG-DAR association table ──────────────────────────────────────────────

deg_dar <- merged_df %>%
  filter(sig_RNA & sig_ATAC) %>%
  dplyr::select(gene_id, label, logFC_RNA, padj_RNA,
                peak_id, logFC_ATAC, padj_ATAC,
                annotation, Direction) %>%
  dplyr::rename(symbol = label) %>%
  arrange(desc(abs(logFC_RNA)))

n_pairs <- nrow(deg_dar)
n_genes <- n_distinct(deg_dar$gene_id)

message(">>> DEG-DAR associations")
message("  Total peak-gene pairs : ", n_pairs)
message("  Unique genes          : ", n_genes)
message("  UP-Opening            : ", sum(deg_dar$Direction == "UP-Opening"))
message("  DOWN-Closing          : ", sum(deg_dar$Direction == "DOWN-Closing"))
message("  UP-Closing            : ", sum(deg_dar$Direction == "UP-Closing"))
message("  DOWN-Opening          : ", sum(deg_dar$Direction == "DOWN-Opening"))

write.csv(deg_dar, PATHS$integration, row.names = FALSE)
message("  -> ", PATHS$integration)

# ── 4. Figure 7A — DEG-DAR association scatter ────────────────────────────────

message(">>> Figure 7A: DEG-DAR scatter")

# Colour and shape per direction — concordant same direction = filled circle,
# discordant = X, matching the visual style of the original figure
dir_colors <- c(
  "UP-Opening"      = "#2c7bb6",
  "DOWN-Closing"    = "#2c7bb6",
  "UP-Closing"      = "#d7191c",
  "DOWN-Opening"    = "#d7191c",
  "Not significant" = "grey80"
)

dir_shapes <- c(
  "UP-Opening"      = 16,
  "DOWN-Closing"    = 16,
  "UP-Closing"      = 16,
  "DOWN-Opening"    = 16,
  "Not significant" = 16
)

dir_sizes <- c(
  "UP-Opening"      = 3.5,
  "DOWN-Closing"    = 3.5,
  "UP-Closing"      = 3.5,
  "DOWN-Opening"    = 3.5,
  "Not significant" = 1.5
)

dir_labels <- c(
  "UP-Opening"      = "Concordant",
  "DOWN-Closing"    = "Concordant",
  "UP-Closing"      = "Opposite direction",
  "DOWN-Opening"    = "Opposite direction",
  "Not significant" = "Not significant"
)

# All DEG-DAR pairs labelled — one arrow per point
label_df <- merged_df %>%
  filter(sig_RNA & sig_ATAC)

# Significant points plotted on top
merged_plot <- merged_df %>%
  arrange(Direction == "Not significant")

lim <- ceiling(max(abs(c(merged_df$logFC_RNA, merged_df$logFC_ATAC)),
                   na.rm = TRUE)) + 0.4

p_scatter <- ggplot(merged_plot, aes(logFC_RNA, logFC_ATAC)) +
  geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.5) +
  geom_vline(xintercept = 0, colour = "grey70", linewidth = 0.5) +
  geom_abline(slope = 1, intercept = 0,
              linetype = "dotted", colour = "grey55", linewidth = 0.5) +
  geom_point(
    aes(color = Direction, shape = Direction, size = Direction),
    alpha = 0.88
  ) +
  scale_color_manual(values = dir_colors, labels = dir_labels, name = NULL) +
  scale_shape_manual(values = dir_shapes, labels = dir_labels, name = NULL) +
  scale_size_manual(values  = dir_sizes,  labels = dir_labels, name = NULL) +
  geom_text_repel(
    data          = label_df,
    aes(label     = label),
    size          = 2.9,
    fontface      = "italic",
    colour        = "grey20",
    box.padding   = 0.35,
    point.padding = 0.2,
    segment.size  = 0.3,
    segment.color = "grey55",
    force         = 2,
    max.overlaps  = Inf,
    show.legend   = FALSE
  ) +
  scale_x_continuous(limits = c(-lim, lim), breaks = -4:4) +
  scale_y_continuous(limits = c(-lim, lim), breaks = -4:4) +
  coord_equal(xlim = c(-lim, lim), ylim = c(-lim, lim)) +
  labs(
    title    = "Concordance of transcriptional and chromatin accessibility changes",
    subtitle = sprintf(
      "%d peak-gene pairs across %d genes (FDR<0.05, |LFC|>%.3f)",
      n_pairs, n_genes, LFC_CUTOFF_SOFT
    ),
    x = expression(RNA~~log[2]~fold~change),
    y = expression(ATAC~~log[2]~fold~change)
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.minor  = element_blank(),
    panel.grid.major  = element_line(colour = "grey93", linewidth = 0.3),
    plot.title        = element_text(face = "bold", size = 11),
    plot.subtitle     = element_text(size = 8.5, colour = "grey40"),
    axis.title        = element_text(size = 10),
    legend.position   = "top",
    legend.spacing.x  = unit(0.3, "cm"),
    legend.text       = element_text(size = 9)
  ) +
  guides(color = guide_legend(override.aes = list(size = 3)))

save_plot(p_scatter, "B04_DEG_DAR_scatter", width = 7, height = 7)

# Focused scatter: sig points only, auto scale, segment arrows on all labels
lim_focused <- max(abs(c(label_df$logFC_RNA, label_df$logFC_ATAC)),
                   na.rm = TRUE) * 1.15

p_scatter_focused <- ggplot(label_df, aes(logFC_RNA, logFC_ATAC)) +
  geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.5) +
  geom_vline(xintercept = 0, colour = "grey70", linewidth = 0.5) +
  geom_abline(slope = 1, intercept = 0,
              linetype = "dotted", colour = "grey55", linewidth = 0.5) +
  geom_point(
    aes(color = Direction, shape = Direction, size = Direction),
    alpha = 0.88
  ) +
  scale_color_manual(values = dir_colors, labels = dir_labels, name = NULL) +
  scale_shape_manual(values = dir_shapes, labels = dir_labels, name = NULL) +
  scale_size_manual(values  = dir_sizes,  labels = dir_labels, name = NULL) +
  geom_text_repel(
    aes(label = label),
    size               = 2.9,
    fontface           = "italic",
    colour             = "grey20",
    box.padding        = 0.4,
    point.padding      = 0.3,
    segment.size       = 0.4,
    segment.color      = "grey40",
    force              = 3,
    min.segment.length = 0,
    max.overlaps       = Inf,
    show.legend        = FALSE
  ) +
  scale_x_continuous(limits = c(-lim_focused, lim_focused)) +
  scale_y_continuous(limits = c(-lim_focused, lim_focused)) +
  coord_equal(
    xlim = c(-lim_focused, lim_focused),
    ylim = c(-lim_focused, lim_focused)
  ) +
  labs(
    title    = "Concordance of transcriptional and chromatin accessibility changes",
    subtitle = sprintf(
      "%d peak-gene pairs across %d genes (FDR<0.05, |LFC|>%.3f)",
      n_pairs, n_genes, LFC_CUTOFF_SOFT
    ),
    x = expression(RNA~~log[2]~fold~change),
    y = expression(ATAC~~log[2]~fold~change)
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(colour = "grey93", linewidth = 0.3),
    plot.title       = element_text(face = "bold", size = 11),
    plot.subtitle    = element_text(size = 8.5, colour = "grey40"),
    axis.title       = element_text(size = 10),
    legend.position  = "top",
    legend.spacing.x = unit(0.3, "cm"),
    legend.text      = element_text(size = 9)
  ) +
  guides(color = guide_legend(override.aes = list(size = 3)))

save_plot(p_scatter_focused, "B04_DEG_DAR_scatter_focused", width = 7, height = 7)


# ── 5. Figure 7B — Dual heatmap ───────────────────────────────────────────────

message(">>> Figure 7B: Dual heatmap")

heatmap_ensg <- unique(deg_dar$gene_id)

# RNA corrected matrix
dds_rna <- readRDS(PATHS$dds_rna)
vsd_rna <- vst(dds_rna, blind = FALSE)
mat_rna_all <- limma::removeBatchEffect(
  assay(vsd_rna),
  covariates = as.matrix(colData(dds_rna)[, paste0("W_", seq_len(K_FACTOR)), drop = FALSE]),
  batch2     = dds_rna$sex,
  design     = model.matrix(~ status, colData(dds_rna))
)
rownames(mat_rna_all) <- sub("\\..*$", "", rownames(mat_rna_all))

ensg_in_rna <- heatmap_ensg[heatmap_ensg %in% rownames(mat_rna_all)]
mat_rna     <- mat_rna_all[ensg_in_rna, , drop = FALSE]

# Row display labels: symbol if available, gene_id otherwise
row_labels <- deg_dar %>%
  distinct(gene_id, symbol) %>%
  filter(gene_id %in% ensg_in_rna) %>%
  group_by(gene_id) %>% slice(1) %>% ungroup()

rownames(mat_rna) <- row_labels$symbol[match(ensg_in_rna, row_labels$gene_id)]

# ATAC corrected matrix
dds_atac <- readRDS(PATHS$dds_atac)
vsd_atac <- vst(dds_atac, blind = FALSE)
mat_atac_all <- limma::removeBatchEffect(
  assay(vsd_atac),
  covariates = as.matrix(colData(dds_atac)[, paste0("W_", seq_len(K_FACTOR)), drop = FALSE]),
  batch2     = dds_atac$sex,
  design     = model.matrix(~ status, colData(dds_atac))
)

# Aggregate peaks per gene: mean VST across all peaks annotated to that ENSG
anno_df <- read.csv(PATHS$atac_annotated) %>%
  dplyr::rename(ensembl_atac = ENSEMBL)

mat_atac_gene <- lapply(ensg_in_rna, function(g) {
  pids <- anno_df %>%
    filter(ensembl_atac == g, peak_id %in% rownames(mat_atac_all)) %>%
    pull(peak_id)
  if (length(pids) == 0) return(NULL)
  colMeans(mat_atac_all[pids, , drop = FALSE])
}) %>% setNames(ensg_in_rna)

keep        <- !sapply(mat_atac_gene, is.null)
ensg_shared <- ensg_in_rna[keep]

if (length(ensg_shared) == 0) stop("No ATAC peaks found for any DEG-DAR gene — check anno_df$ensembl_atac vs ensg_in_rna")

mat_atac <- do.call(rbind, mat_atac_gene[keep])
shared_labels <- row_labels$symbol[match(ensg_shared, row_labels$gene_id)]
rownames(mat_rna)  <- shared_labels[match(ensg_in_rna,  ensg_shared)]
rownames(mat_atac) <- shared_labels
mat_rna <- mat_rna[match(ensg_shared, ensg_in_rna), , drop = FALSE]

message("  Genes in dual heatmap: ", length(ensg_shared))

# Shared row order from RNA clustering
row_order    <- hclust(dist(t(scale(t(mat_rna)))), method = "ward.D2")$order
mat_rna_ord  <- mat_rna[row_order,  , drop = FALSE]
mat_atac_ord <- mat_atac[row_order, , drop = FALSE]

# Column order: status then sex
meta_rna_df  <- as.data.frame(colData(dds_rna))[,  c("status", "sex")]
meta_atac_df <- as.data.frame(colData(dds_atac))[, c("status", "sex")]
col_ord_rna  <- order(meta_rna_df$status,  meta_rna_df$sex)
col_ord_atac <- order(meta_atac_df$status, meta_atac_df$sex)
mat_rna_ord  <- mat_rna_ord[,  col_ord_rna,  drop = FALSE]
mat_atac_ord <- mat_atac_ord[, col_ord_atac, drop = FALSE]
anno_col_rna  <- meta_rna_df[col_ord_rna,  , drop = FALSE]
anno_col_atac <- meta_atac_df[col_ord_atac, , drop = FALSE]
anno_colors   <- list(status = COLORS_STATUS, sex = COLORS_SEX)

fsize_row <- if (length(ensg_shared) > 40) 5 else if (length(ensg_shared) > 20) 7 else 9

ph_rna <- pheatmap(
  mat_rna_ord,
  scale             = "row",
  color             = colorRampPalette(c("navy", "white", "firebrick3"))(100),
  breaks            = seq(-3, 3, length.out = 101),
  cluster_rows      = FALSE,
  cluster_cols      = FALSE,
  annotation_col    = anno_col_rna,
  annotation_colors = anno_colors,
  show_colnames     = FALSE,
  show_rownames     = TRUE,
  fontsize_row      = fsize_row,
  main              = "RNA-seq (VST, z-scored)",
  legend            = TRUE,
  silent            = TRUE
)

ph_atac <- pheatmap(
  mat_atac_ord,
  scale             = "row",
  color             = colorRampPalette(c("navy", "white", "firebrick3"))(100),
  breaks            = seq(-3, 3, length.out = 101),
  cluster_rows      = FALSE,
  cluster_cols      = FALSE,
  annotation_col    = anno_col_atac,
  annotation_colors = anno_colors,
  show_colnames     = FALSE,
  show_rownames     = FALSE,
  main              = "ATAC-seq (mean VST per gene, z-scored)",
  legend            = FALSE,
  silent            = TRUE
)

pdf(
  file.path(PATHS$plots, "B04_DEG_DAR_heatmap.pdf"),
  width  = 14,
  height = max(6, length(ensg_shared) * 0.25 + 3)
)
grid.arrange(ph_rna$gtable, ph_atac$gtable,
             ncol = 2, widths = c(1.4, 1))
dev.off()
message("  -> B04_DEG_DAR_heatmap.pdf")

# ── 6. GO enrichment on DEG-DAR genes ────────────────────────────────────────

message(">>> GO enrichment: DEG-DAR genes")

go_symbols <- deg_dar %>%
  filter(!is.na(symbol), symbol != "", !grepl("^ENSG", symbol)) %>%
  pull(symbol) %>% unique()

gene_ids <- bitr(go_symbols, fromType = "SYMBOL", toType = "ENTREZID",
                 OrgDb = org.Hs.eg.db, drop = TRUE)

ego <- enrichGO(
  gene          = gene_ids$ENTREZID,
  OrgDb         = org.Hs.eg.db,
  ont           = GO_ONT,
  pAdjustMethod = "BH",
  pvalueCutoff  = P_CUTOFF,
  readable      = TRUE
)

if (!is.null(ego) && nrow(ego) > 0) {
  ego_s <- clusterProfiler::simplify(ego, cutoff = GO_SIMPLIFY,
                                     by = "p.adjust", select_fun = min)
  write.csv(
    as.data.frame(ego_s),
    file.path(PATHS$tables, "B04_DEG_DAR_GO_BP.csv"),
    row.names = FALSE
  )
  p_go <- dotplot(ego_s, showCategory = 20) +
    labs(
      title    = "GO Biological Process — DEG-DAR genes",
      subtitle = sprintf("Input: %d genes | p.adj < %s (BH)", length(go_symbols), P_CUTOFF)
    ) +
    theme(plot.title = element_text(size = 11))
  save_plot(p_go, "B04_DEG_DAR_GO_dotplot", width = 9, height = 10)
  message("  -> B04_DEG_DAR_GO_BP.csv")
} else {
  message("  No significant GO terms")
}

# ── 7. Summary ────────────────────────────────────────────────────────────────

message("\n>>> Summary")
message("  DEG-DAR pairs   : ", n_pairs)
message("  Unique genes    : ", n_genes)
message("  Heatmap genes   : ", length(ensg_shared))
message("  UP-Opening      : ", sum(deg_dar$Direction == "UP-Opening"))
message("  DOWN-Closing    : ", sum(deg_dar$Direction == "DOWN-Closing"))
message("  UP-Closing      : ", sum(deg_dar$Direction == "UP-Closing"))
message("  DOWN-Opening    : ", sum(deg_dar$Direction == "DOWN-Opening"))

save_session("B04")
message("\n✓ B04 complete.\n")