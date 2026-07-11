# B06_figures.R
# Publication figure generation

source("scripts/00_config.R")
init_dirs()

suppressPackageStartupMessages({
  library(tidyverse); library(ggrepel)
  library(DESeq2); library(limma); library(pheatmap)
})

# Base theme for publication figures
BASE_THEME <- theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(colour = "grey93", linewidth = 0.3),
        plot.title = element_text(face = "bold", size = 11),
        plot.subtitle = element_text(size = 8.5, colour = "grey40"),
        axis.title = element_text(size = 10),
        legend.text = element_text(size = 9))

STATUS_COLORS <- c(control = "#4575b4", KS_I = "#d73027")
STATUS_LABELS <- c(control = "Control", KS_I = "KS type I")

# PCA calculation helper
calc_pca <- function(mat, meta) {
  pca <- prcomp(t(mat))
  df <- as.data.frame(pca$x) %>%
    rownames_to_column("sample_id") %>%
    left_join(meta, by = "sample_id")
  var <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)
  list(df = df, var = var)
}

# Add display columns for plotting
add_display_cols <- function(df) {
  df %>% mutate(
    grp = case_when(
      as.character(status) == "control" ~ "Control",
      as.character(status) == "KS_I" ~ "KS type I",
      TRUE ~ as.character(status)
    ),
    grp = factor(grp, levels = c("Control", "KS type I")),
    bat = factor(paste0("Batch ", as.character(batch)),
                 levels = c("Batch 84", "Batch 185"))
  )
}

# Load RNA-seq data
dds_rna <- readRDS(PATHS$dds_rna)
vsd_rna <- vst(dds_rna, blind = FALSE)
meta_rna <- as.data.frame(colData(dds_rna)) %>%
  rownames_to_column("sample_id") %>%
  mutate(batch = as.character(batch))

# Raw and corrected expression matrices
mat_raw <- assay(vsd_rna)
mat_cor <- limma::removeBatchEffect(
  mat_raw,
  covariates = as.matrix(colData(dds_rna)[, paste0("W_", seq_len(K_FACTOR)), drop = FALSE]),
  batch2 = dds_rna$sex,
  design = model.matrix(~ status, colData(dds_rna))
)

pca_raw <- calc_pca(mat_raw, meta_rna)
pca_cor <- calc_pca(mat_cor, meta_rna)
pca_raw$df <- add_display_cols(pca_raw$df)
pca_cor$df <- add_display_cols(pca_cor$df)

# PCA: uncorrected
p_pca_raw <- ggplot(pca_raw$df, aes(x = PC1, y = PC2, colour = grp, shape = bat)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey75") +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey75") +
  geom_point(size = 3.5, alpha = 0.85) +
  scale_colour_manual(values = c("Control" = "#4575b4", "KS type I" = "#d73027"),
                      name = "Disease status") +
  scale_shape_manual(values = c("Batch 84" = 16L, "Batch 185" = 17L),
                     name = "Sequencing batch") +
  labs(title = "PCA: RNA-seq counts (uncorrected)",
       x = paste0("PC1: ", pca_raw$var[1], "% variance"),
       y = paste0("PC2: ", pca_raw$var[2], "% variance")) +
  BASE_THEME
save_plot(p_pca_raw, "PCA_raw", width = 6.5, height = 5)

# PCA: RUV-corrected
p_pca_cor <- ggplot(pca_cor$df, aes(x = PC1, y = PC2, colour = grp, shape = bat)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey75") +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey75") +
  geom_point(size = 3.5, alpha = 0.85) +
  scale_colour_manual(values = c("Control" = "#4575b4", "KS type I" = "#d73027"),
                      name = "Disease status") +
  scale_shape_manual(values = c("Batch 84" = 16L, "Batch 185" = 17L),
                     name = "Sequencing batch") +
  labs(title = "PCA: RNA-seq counts (RUV-corrected)",
       x = paste0("PC1: ", pca_cor$var[1], "% variance"),
       y = paste0("PC2: ", pca_cor$var[2], "% variance")) +
  BASE_THEME
save_plot(p_pca_cor, "PCA_corrected", width = 6.5, height = 5)

# Relative Log Expression (RLE)
rle_df <- mat_cor %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  pivot_longer(-gene, names_to = "sample_id", values_to = "expression") %>%
  group_by(gene) %>%
  mutate(rle = expression - median(expression)) %>%
  ungroup() %>%
  left_join(meta_rna, by = "sample_id") %>%
  mutate(sample_label = sub("_b(84|185)$", "", sample_id))

sample_order <- rle_df %>%
  distinct(sample_label, status) %>%
  arrange(status, sample_label) %>%
  pull(sample_label)

rle_df <- rle_df %>%
  mutate(sample_label = factor(sample_label, levels = sample_order))

p_rle <- ggplot(rle_df, aes(x = sample_label, y = rle, fill = status)) +
  geom_boxplot(outlier.size = 0.4, outlier.alpha = 0.25,
               alpha = 0.8, colour = "grey25", linewidth = 0.3) +
  geom_hline(yintercept = 0, colour = "firebrick", linetype = "dashed", linewidth = 0.5) +
  scale_fill_manual(values = STATUS_COLORS, labels = STATUS_LABELS, name = "Disease status") +
  labs(title = "Relative Log Expression",
       x = "Sample", y = "Deviation from gene median") +
  BASE_THEME +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))
save_plot(p_rle, "RLE", width = 8, height = 4)

# Volcano plot
volcano_df <- read.csv(PATHS$rna_de) %>%
  filter(!is.na(padj)) %>%
  mutate(
    significance = case_when(
      padj < P_CUTOFF & log2FoldChange > LFC_CUTOFF ~ "Upregulated",
      padj < P_CUTOFF & log2FoldChange < -LFC_CUTOFF ~ "Downregulated",
      TRUE ~ "Not significant"
    ),
    is_artifact = (symbol == "NPIPA8"),
    log10_padj = pmin(-log10(padj), 50)
  )

n_up <- sum(volcano_df$significance == "Upregulated")
n_down <- sum(volcano_df$significance == "Downregulated")

label_genes <- volcano_df %>%
  filter(!is_artifact, significance != "Not significant") %>%
  slice_min(padj, n = 10) %>%
  pull(symbol)

p_volcano <- ggplot(volcano_df, aes(x = log2FoldChange, y = log10_padj)) +
  geom_vline(xintercept = c(-LFC_CUTOFF, LFC_CUTOFF),
             linetype = "dashed", colour = "grey55", linewidth = 0.4) +
  geom_hline(yintercept = -log10(P_CUTOFF),
             linetype = "dashed", colour = "grey55", linewidth = 0.4) +
  geom_point(data = filter(volcano_df, significance == "Not significant"),
             colour = "grey85", size = 1.2, alpha = 0.5) +
  geom_point(data = filter(volcano_df, significance != "Not significant", !is_artifact),
             aes(colour = significance), size = 1.8, alpha = 0.75) +
  geom_point(data = filter(volcano_df, is_artifact),
             shape = 2, size = 3, colour = "black", stroke = 1) +
  geom_text_repel(data = filter(volcano_df, symbol %in% label_genes),
                  aes(label = symbol), size = 3.2, fontface = "italic",
                  box.padding = 0.5, segment.color = "grey45", max.overlaps = 20) +
  scale_colour_manual(values = c("Upregulated" = "#d73027", "Downregulated" = "#4575b4"),
                      name = NULL) +
  labs(title = "Differential gene expression",
       subtitle = sprintf("%d upregulated, %d downregulated", n_up, n_down),
       x = expression(log[2]~fold~change),
       y = expression(-log[10]~FDR)) +
  BASE_THEME + theme(legend.position = "top")
save_plot(p_volcano, "volcano", width = 7, height = 6)

# ── Fig 3B: DEG heatmap ───────────────────────────────────────────────────────
message("> Fig 3B: DEG heatmap (top 100, |LFC| > 1)...")

top_degs <- read.csv(PATHS$rna_de) %>%
  filter(!is.na(padj),
         padj < P_CUTOFF,
         abs(log2FoldChange) > LFC_CUTOFF,
         symbol != "NPIPA8",
         symbol != "") %>%
  slice_min(padj, n = 100, with_ties = FALSE) %>%
  pull(symbol)

gene_map <- read.table(PATHS$gene_map, header = TRUE, sep = "\t",
                       row.names = 1, check.names = FALSE) %>%
  rownames_to_column("gene_id")   # Ensembl ID is the rowname in B01 output

ens_ids <- gene_map %>%
  filter(gene_name %in% top_degs) %>%
  dplyr::select(gene_id, gene_name) %>%
  distinct(gene_name, .keep_all = TRUE)

heatmap_mat <- mat_cor[
  intersect(ens_ids$gene_id, rownames(mat_cor)), , drop = FALSE
]
rownames(heatmap_mat) <- ens_ids$gene_name[
  match(rownames(heatmap_mat), ens_ids$gene_id)
]

# Z-score rows; cap at ±3 to prevent outliers from washing out the palette
heatmap_z <- t(scale(t(heatmap_mat)))
heatmap_z <- pmax(pmin(heatmap_z, 3), -3)

anno_col <- meta_rna %>%
  dplyr::select(sample_id, status) %>%
  column_to_rownames("sample_id")

anno_colors <- list(
  status = setNames(STATUS_COLORS, names(STATUS_COLORS))
)

# Column order: controls then KS, alphabetical within group
col_order <- meta_rna %>%
  arrange(status, sample_id) %>%
  pull(sample_id)
heatmap_z <- heatmap_z[, intersect(col_order, colnames(heatmap_z)), drop = FALSE]

pheatmap::pheatmap(
  heatmap_z,
  annotation_col    = anno_col,
  annotation_colors = anno_colors,
  color             = colorRampPalette(c("#4575b4", "white", "#d73027"))(100),
  cluster_cols      = FALSE,
  cluster_rows      = TRUE,
  show_colnames     = FALSE,
  fontsize_row      = 6,
  border_color      = NA,
  main              = sprintf("Top %d DEGs — KS type I vs control (z-score, capped ±3)",
                              nrow(heatmap_z)),
  filename          = file.path(PATHS$plots, "B06_Fig3B_DEG_heatmap.pdf"),
  width             = 7,
  height            = 14
)
message("  -> ", file.path(PATHS$plots, "B06_Fig3B_DEG_heatmap.pdf"))











# ── GO enrichment plots (Reviewer 2) ───────────────────────────────
TABLES <- PATHS$tables

# Explicit column types: a GO table can legitimately be header-only with
# zero rows (nothing significant under the current universe/threshold).
# read.csv() can't infer types from empty data and defaults such columns
# to logical, which breaks bind_rows() against a non-empty table of the
# same shape — so pin the enrichGO() output schema explicitly.
GO_COL_TYPES <- c(ID = "character", Description = "character",
                  GeneRatio = "character", BgRatio = "character",
                  pvalue = "numeric", p.adjust = "numeric",
                  qvalue = "numeric", geneID = "character",
                  Count = "integer")
read_go_csv <- function(path) read.csv(path, colClasses = GO_COL_TYPES)

parse_ratio <- function(x) {
  vapply(x, function(r) {
    v <- as.numeric(strsplit(r, "/")[[1]])
    v[1] / v[2]
  }, numeric(1))
}

go_dotplot_revised <- function(go_df, title, subtitle, n = 15) {
  df <- go_df %>%
    filter(!is.na(p.adjust)) %>%
    arrange(p.adjust) %>%
    slice_head(n = n) %>%
    mutate(
      enrichment_ratio = parse_ratio(GeneRatio) / parse_ratio(BgRatio),
      label = factor(
        str_wrap(Description, 40),
        levels = str_wrap(Description[order(-enrichment_ratio)], 40)
      )
    )
  
  ggplot(df, aes(x = enrichment_ratio, y = label, colour = -log10(p.adjust))) +
    geom_point(size = 8, alpha = 0.9) +
    geom_text(aes(label = Count), colour = "white", size = 2.8, fontface = "bold") +
    scale_colour_gradientn(
      colours = c("#2166ac", "#4393c3", "#d6604d", "#b2182b"),
      name    = expression(-log[10](FDR)),
      labels  = scales::label_number(accuracy = 0.1)
    ) +
    labs(title = title, subtitle = subtitle,
         x = "Enrichment ratio", y = NULL) +
    theme_bw(base_size = 11) +
    theme(
      panel.grid.minor  = element_blank(),
      panel.grid.major  = element_line(colour = "grey93", linewidth = 0.3),
      plot.title        = element_text(face = "bold", size = 11),
      plot.subtitle     = element_text(size = 8.5, colour = "grey40"),
      axis.text.y       = element_text(size = 9)
    )
}

# Fig 4B — RNA GO
go_rna_strict <- read_go_csv(file.path(TABLES, "B02_RNA_GO_BP_strict.csv"))

p_rna_go_revised <- go_dotplot_revised(
  go_rna_strict,
  title    = "GO Biological Process enrichment",
  subtitle = "Top 15 terms | FDR < 0.05, |LFC| > 1"
)
save_plot(p_rna_go_revised, "B06_Fig4B_RNA_GO_revised", width = 7, height = 7)

# Fig 6C/D — ATAC GO
go_open  <- read_go_csv(file.path(TABLES, "B03_ATAC_GO_opening.csv"))
go_close <- read_go_csv(file.path(TABLES, "B03_ATAC_GO_closing.csv"))

go_atac_combined <- bind_rows(
  go_open  %>% mutate(direction = "Opening"),
  go_close %>% mutate(direction = "Closing")
) %>%
  filter(!is.na(p.adjust)) %>%
  group_by(direction) %>%
  arrange(p.adjust) %>%
  slice_head(n = 15) %>%
  ungroup() %>%
  mutate(
    direction        = factor(direction, levels = c("Opening", "Closing")),
    enrichment_ratio = parse_ratio(GeneRatio) / parse_ratio(BgRatio),
    label            = str_wrap(Description, 35)
  )

p_atac_go_revised <- ggplot(go_atac_combined,
                            aes(x = enrichment_ratio,
                                y = reorder(label, enrichment_ratio),
                                colour = -log10(p.adjust))) +
  geom_point(size = 8, alpha = 0.9) +
  geom_text(aes(label = Count), colour = "white", size = 2.8, fontface = "bold") +
  scale_colour_gradientn(
    colours = c("#2166ac", "#4393c3", "#d6604d", "#b2182b"),
    name    = expression(-log[10](FDR)),
    labels  = scales::label_number(accuracy = 0.1)
  ) +
  facet_wrap(~ direction, scales = "free", ncol = 2) +
  labs(
    title    = "GO Biological Process enrichment: ATAC-seq DA peaks",
    subtitle = "Top 15 terms per direction | FDR < 0.05, |LFC| > 1",
    x = "Enrichment ratio", y = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.minor  = element_blank(),
    panel.grid.major  = element_line(colour = "grey93", linewidth = 0.3),
    strip.background  = element_rect(fill = "grey95", colour = "grey80"),
    strip.text        = element_text(face = "bold", size = 10),
    plot.title        = element_text(face = "bold", size = 11),
    plot.subtitle     = element_text(size = 8.5, colour = "grey40"),
    axis.text.y       = element_text(size = 8.5)
  )
save_plot(p_atac_go_revised, "B06_Fig6CD_ATAC_GO_revised", width = 13, height = 7)

# Fig 7B — Integration GO
go_int <- read_go_csv(file.path(TABLES, "B04_DEG_DAR_GO_BP.csv"))

p_int_go_revised <- go_dotplot_revised(
  go_int,
  title    = "GO Biological Process enrichment: DEG-DAR genes",
  subtitle = sprintf("All %d significant terms | FDR < 0.05", nrow(go_int %>% filter(!is.na(p.adjust))))
)
save_plot(p_int_go_revised, "B06_Fig7B_integration_GO_revised", width = 7, height = 5)

save_session("B06")