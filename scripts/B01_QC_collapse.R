# B01_QC_collapse.R
# Quality control and technical replicate collapsing

source("00_config.R")
init_dirs()

suppressPackageStartupMessages({
  library(tidyverse)
  library(janitor)
  library(DESeq2)
  library(limma)
  library(pheatmap)
  library(ggrepel)
  library(ggpubr)
})

# Load raw counts and metadata
raw_full <- read.table(PATHS$counts_raw,
                       header = TRUE, row.names = 1, check.names = FALSE,
                       sep = "\t")

if (!"gene_name" %in% colnames(raw_full)) stop("Column 'gene_name' not found.")
gene_map   <- raw_full[, "gene_name", drop = FALSE]
raw_counts <- raw_full %>% dplyr::select(-gene_name)

ss_full <- read.csv(PATHS$sample_sheet, stringsAsFactors = FALSE) %>%
  janitor::clean_names() %>%
  dplyr::mutate(
    status = factor(status, levels = c("control", "KS_I")),
    sex    = factor(sex),
    batch  = factor(batch)
  )

# Align sample sheet with count matrix
common_libs <- intersect(ss_full$analysis_id, colnames(raw_counts))
ss_full     <- ss_full %>% dplyr::filter(analysis_id %in% common_libs)
raw_counts  <- raw_counts[, ss_full$analysis_id]

message(sprintf("  Loaded: %d genes x %d libraries across %d GDBs",
                nrow(raw_counts),
                ncol(raw_counts),
                length(unique(ss_full$gdb))))

# QC on all libraries before collapsing
rownames(ss_full) <- ss_full$analysis_id

dds_raw <- DESeqDataSetFromMatrix(
  countData = round(as.matrix(raw_counts)),
  colData   = ss_full,
  design    = ~ batch + status
)
dds_raw <- dds_raw[rowSums(counts(dds_raw) >= 10) >= 3, ]
vsd_raw <- vst(dds_raw, blind = TRUE)

colnames(vsd_raw) <- paste0(colData(vsd_raw)$batch, "_", colData(vsd_raw)$gdb)

# PCA: uncorrected
pcaData    <- DESeq2::plotPCA(vsd_raw, intgroup = c("status", "batch"),
                               returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
pcaData$label <- colnames(vsd_raw)

p_pca_raw <- ggplot(pcaData, aes(PC1, PC2, color = status, shape = batch)) +
  geom_point(size = 3) +
  geom_text_repel(aes(label = label), size = 2.5, max.overlaps = Inf,
                  box.padding = 0.4, show.legend = FALSE) +
  scale_color_manual(values = COLORS_STATUS) +
  scale_shape_manual(values = c("84" = 16, "185" = 17)) +
  labs(title = "PCA all libraries (pre-collapse, blind VST)",
       x = paste0("PC1: ", percentVar[1], "% variance"),
       y = paste0("PC2: ", percentVar[2], "% variance")) +
  theme_bw() + theme(aspect.ratio = 1)

save_plot(p_pca_raw, "B01_QC_PCA_raw")

# PCA: batch and sex corrected (visualization only)
design_protect <- model.matrix(~ status, data = as.data.frame(colData(vsd_raw)))
mat_corr       <- limma::removeBatchEffect(
  assay(vsd_raw),
  batch  = vsd_raw$batch,
  batch2 = vsd_raw$sex,
  design = design_protect
)
vsd_corr        <- vsd_raw
assay(vsd_corr) <- mat_corr

pcaData2    <- DESeq2::plotPCA(vsd_corr, intgroup = c("status", "batch"),
                                returnData = TRUE)
percentVar2 <- round(100 * attr(pcaData2, "percentVar"))
pcaData2$label <- colnames(vsd_corr)

p_pca_corr <- ggplot(pcaData2, aes(PC1, PC2, color = status, shape = batch)) +
  geom_point(size = 3) +
  geom_text_repel(aes(label = label), size = 2.5, max.overlaps = Inf,
                  box.padding = 0.4, show.legend = FALSE) +
  scale_color_manual(values = COLORS_STATUS) +
  scale_shape_manual(values = c("84" = 16, "185" = 17)) +
  labs(title = "PCA all libraries (batch + sex corrected, visualisation only)",
       x = paste0("PC1: ", percentVar2[1], "% variance"),
       y = paste0("PC2: ", percentVar2[2], "% variance")) +
  theme_bw() + theme(aspect.ratio = 1)

save_plot(p_pca_corr, "B01_QC_PCA_corrected")

# Sample distance heatmap
sampleDistMatrix <- as.matrix(dist(t(mat_corr)))
anno_col <- as.data.frame(colData(vsd_corr)[, c("status", "batch", "sex")])

save_heatmap(
  mat         = sampleDistMatrix,
  anno_col    = anno_col,
  anno_colors = list(status = COLORS_STATUS,
                     sex    = COLORS_SEX,
                     batch  = COLORS_BATCH),
  title       = "Sample distances (batch + sex corrected)",
  name        = "B01_QC_sample_distances",
  show_rownames = TRUE,
  fontsize_row  = 7
)

# Library size distribution
lib_sizes <- data.frame(
  library   = colnames(raw_counts),
  total     = colSums(raw_counts) / 1e6,
  stringsAsFactors = FALSE
) %>%
  left_join(ss_full %>% dplyr::select(analysis_id, status, batch, sex),
            by = c("library" = "analysis_id")) %>%
  arrange(total)

p_libsize <- ggplot(lib_sizes,
                    aes(x = reorder(library, total), y = total, fill = status)) +
  geom_col() +
  geom_hline(yintercept = 20, linetype = "dashed", color = "grey40") +
  scale_fill_manual(values = COLORS_STATUS) +
  coord_flip() +
  labs(title = "Library sizes (raw mapped reads)",
       x = NULL, y = "Total counts (millions)", fill = "Status") +
  theme_bw(base_size = 10)

save_plot(p_libsize, "B01_QC_library_sizes", width = 8, height = 6)

# Per-gene count distribution
mean_counts <- rowMeans(as.matrix(raw_counts))
p_countdist <- ggplot(data.frame(mean = log10(mean_counts + 1)),
                      aes(x = mean)) +
  geom_histogram(bins = 80, fill = "steelblue", color = "white", linewidth = 0.2) +
  geom_vline(xintercept = log10(10), linetype = "dashed", color = "firebrick") +
  annotate("text", x = log10(10) + 0.05, y = Inf,
           label = "filter threshold (>=10)", vjust = 2, hjust = 0,
           size = 3, color = "firebrick") +
  labs(title = "Gene count distribution (all libraries)",
       x = "log10(mean raw count + 1)", y = "Number of genes") +
  theme_bw()

save_plot(p_countdist, "B01_QC_count_distribution", width = 7, height = 5)

# Technical replicate consistency (same GDB, different batch)
df_meta    <- as.data.frame(colData(vsd_raw)) %>% rownames_to_column("lib_id")
rep_groups <- df_meta %>%
  group_by(gdb) %>%
  filter(n() > 1) %>%
  summarise(libs = list(lib_id), .groups = "drop")

if (nrow(rep_groups) > 0) {
  plot_data <- map_dfr(seq_len(nrow(rep_groups)), function(i) {
    libs <- unlist(rep_groups$libs[[i]])
    if (length(libs) < 2) return(NULL)
    vals <- assay(vsd_raw)[, libs[1:2]]
    data.frame(gdb = rep_groups$gdb[i], Rep1 = vals[, 1], Rep2 = vals[, 2])
  })

  p_reps <- ggplot(plot_data, aes(Rep1, Rep2)) +
    geom_point(alpha = 0.15, size = 0.6, color = "navy") +
    geom_abline(slope = 1, intercept = 0,
                color = "firebrick", linetype = "dashed") +
    stat_cor(method = "pearson", size = 2.8) +
    facet_wrap(~ gdb, ncol = 3) +
    theme_bw() +
    labs(title = "Technical replicate consistency (VST)",
         subtitle = "Each panel: same GDB sequenced in batch 84 and batch 185",
         x = "Batch 84 (VST)", y = "Batch 185 (VST)")

  save_plot(p_reps, "B01_QC_replicate_scatter", width = 10, height = 5)
}

# Detected genes per library (count >= 10)
detected <- apply(as.matrix(raw_counts), 2, function(x) sum(x >= 10))
det_df   <- data.frame(
  library  = names(detected),
  detected = detected,
  stringsAsFactors = FALSE
) %>%
  left_join(ss_full %>% dplyr::select(analysis_id, status, batch),
            by = c("library" = "analysis_id"))

p_detected <- ggplot(det_df,
                     aes(x = reorder(library, detected), y = detected / 1e3,
                         fill = status)) +
  geom_col() +
  scale_fill_manual(values = COLORS_STATUS) +
  coord_flip() +
  labs(title = "Detected genes per library (count >= 10)",
       x = NULL, y = "Genes detected (thousands)", fill = "Status") +
  theme_bw(base_size = 10)

save_plot(p_detected, "B01_QC_detected_genes", width = 8, height = 6)

# Collapse technical replicates by summing counts per GDB
counts_t <- as.data.frame(t(raw_counts)) %>%
  rownames_to_column("analysis_id") %>%
  left_join(dplyr::select(ss_full, analysis_id, gdb), by = "analysis_id")

collapsed_counts <- counts_t %>%
  dplyr::select(-analysis_id) %>%
  group_by(gdb) %>%
  summarise(across(where(is.numeric), sum), .groups = "drop") %>%
  column_to_rownames("gdb") %>%
  t() %>%
  as.data.frame()

# Collapse metadata: take consensus values per GDB
collapsed_meta <- ss_full %>%
  group_by(gdb) %>%
  summarise(across(everything(), function(x) {
    u <- unique(na.omit(as.character(x)))
    if (length(u) == 0) return(NA_character_)
    if (length(u) == 1) return(u)
    paste(sort(u), collapse = ";")
  }), .groups = "drop") %>%
  dplyr::rename(GDB = gdb)

meta_rna    <- load_cohort("RNA")
RNA_SAMPLES <- rownames(meta_rna)

collapsed_counts <- collapsed_counts[, colnames(collapsed_counts) %in% RNA_SAMPLES]
collapsed_meta   <- collapsed_meta %>% dplyr::filter(GDB %in% RNA_SAMPLES)

message(sprintf("  Collapsed: %d genes x %d biological samples",
                nrow(collapsed_counts), ncol(collapsed_counts)))
message("  GDBs: ", paste(sort(colnames(collapsed_counts)), collapse = ", "))

# QC on collapsed data
ss_collapsed <- meta_rna
rownames(ss_collapsed) <- ss_collapsed$gdb

dds_collapsed <- DESeqDataSetFromMatrix(
  countData = round(as.matrix(collapsed_counts[, rownames(ss_collapsed)])),
  colData   = ss_collapsed,
  design    = ~ batch + status
)
dds_collapsed <- dds_collapsed[rowSums(counts(dds_collapsed) >= 10) >= 3, ]
vsd_collapsed <- vst(dds_collapsed, blind = TRUE)

pcaC    <- DESeq2::plotPCA(vsd_collapsed, intgroup = c("status", "batch"),
                            returnData = TRUE)
pctC    <- round(100 * attr(pcaC, "percentVar"))
pcaC$label <- rownames(pcaC)

p_pca_collapsed <- ggplot(pcaC, aes(PC1, PC2, color = status, shape = batch)) +
  geom_point(size = 3.5) +
  geom_text_repel(aes(label = label), size = 2.8, max.overlaps = Inf,
                  box.padding = 0.4, show.legend = FALSE) +
  scale_color_manual(values = COLORS_STATUS) +
  scale_shape_manual(values = c("84" = 16, "185" = 17)) +
  labs(title = "PCA after replicate collapse (blind VST)",
       x = paste0("PC1: ", pctC[1], "% variance"),
       y = paste0("PC2: ", pctC[2], "% variance")) +
  theme_bw() + theme(aspect.ratio = 1)

save_plot(p_pca_collapsed, "B01_QC_PCA_collapsed")

# Save outputs
write.table(collapsed_counts, PATHS$counts_collapsed,
            sep = "\t", quote = FALSE, col.names = NA)
write.csv(collapsed_meta, PATHS$meta_collapsed, row.names = FALSE)
write.table(gene_map, PATHS$gene_map,
            sep = "\t", quote = FALSE, col.names = NA)

save_session("B01")
