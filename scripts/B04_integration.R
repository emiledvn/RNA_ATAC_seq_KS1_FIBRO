# B04_integration.R
# Integration of RNA-seq and ATAC-seq results

source("00_config.R")
init_dirs()

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggrepel)
})

# Load RNA differential expression results
res_rna <- read.csv(PATHS$rna_de) %>%
  dplyr::select(symbol, log2FoldChange, padj) %>%
  dplyr::rename(logFC_RNA = log2FoldChange, padj_RNA = padj) %>%
  filter(!is.na(symbol), symbol != "") %>%
  group_by(symbol) %>% slice_min(padj_RNA, n = 1, with_ties = FALSE) %>% ungroup()

# Load ATAC differential accessibility results (annotated peaks)
res_atac <- read.csv(PATHS$atac_annotated) %>%
  dplyr::rename(symbol = SYMBOL) %>%
  dplyr::select(symbol, log2FoldChange, padj, annotation) %>%
  dplyr::rename(logFC_ATAC = log2FoldChange, padj_ATAC = padj) %>%
  filter(!is.na(symbol)) %>%
  group_by(symbol) %>% slice_min(padj_ATAC, n = 1, with_ties = FALSE) %>% ungroup()

# Merge RNA and ATAC results
merged_df <- inner_join(res_rna, res_atac, by = "symbol") %>%
  mutate(
    sig_RNA  = padj_RNA  < P_CUTOFF & abs(logFC_RNA)  > LFC_CUTOFF_SOFT,
    sig_ATAC = padj_ATAC < P_CUTOFF & abs(logFC_ATAC) > LFC_CUTOFF_SOFT,
    Category = case_when(
      sig_RNA & sig_ATAC & sign(logFC_RNA) == sign(logFC_ATAC) & logFC_RNA > 0 ~ "Concordant_UP",
      sig_RNA & sig_ATAC & sign(logFC_RNA) == sign(logFC_ATAC) & logFC_RNA < 0 ~ "Concordant_DOWN",
      sig_RNA & sig_ATAC & sign(logFC_RNA) != sign(logFC_ATAC)                 ~ "Discordant",
      TRUE ~ "Not_significant"
    )
  )

# Correlation between RNA and ATAC fold changes
cor_res    <- cor.test(merged_df$logFC_RNA, merged_df$logFC_ATAC, method = "pearson")
concordant <- merged_df %>%
  filter(Category %in% c("Concordant_UP", "Concordant_DOWN")) %>%
  arrange(desc(abs(logFC_RNA)))

write.csv(concordant, PATHS$integration, row.names = FALSE)

# Export all significant overlapping genes
overlap_soft <- merged_df %>%
  filter(sig_RNA & sig_ATAC) %>%
  mutate(
    rna_lfc   = logFC_RNA, atac_lfc  = logFC_ATAC,
    rna_padj  = padj_RNA,  atac_padj = padj_ATAC,
    concordant = Category %in% c("Concordant_UP", "Concordant_DOWN")
  ) %>%
  dplyr::select(symbol, rna_lfc, rna_padj, atac_lfc, atac_padj,
                annotation, concordant, Category)

write.csv(overlap_soft,
          file.path(PATHS$tables, "B04_RNA_ATAC_overlap_soft.csv"),
          row.names = FALSE)

message("  Genes in merged dataset : ", nrow(merged_df))
message("  Concordant UP           : ", sum(merged_df$Category == "Concordant_UP"))
message("  Concordant DOWN         : ", sum(merged_df$Category == "Concordant_DOWN"))
message("  Discordant              : ", sum(merged_df$Category == "Discordant"))
message("  Pearson R (all genes)   : ", round(cor_res$estimate, 3))
message("  Pearson p-value         : ", formatC(cor_res$p.value, format = "e", digits = 2))

# Candidate genes of interest (Wnt/PCP pathway)
wnt_pcp <- c("VANGL2","WNT5A","ROR2","CELSR1","PRICKLE1",
             "FZD6","DVL1","DVL3","RORB","ANKRD6","DAAM1")

print(merged_df %>%
  filter(symbol %in% wnt_pcp) %>%
  dplyr::select(symbol, logFC_RNA, padj_RNA, logFC_ATAC, padj_ATAC, Category), n = Inf)

# Quadrant plot: RNA vs ATAC log2 fold changes
message(">>> Quadrant plot...")

top_labels   <- head(concordant$symbol, 15)
story_labels <- intersect(wnt_pcp, merged_df$symbol)
all_labels   <- unique(c(top_labels, story_labels))

cat_colors <- c(Concordant_UP = "firebrick", Concordant_DOWN = "navy",
                Discordant = "darkorange", Not_significant = "grey88")

p_quad <- ggplot(merged_df, aes(logFC_ATAC, logFC_RNA)) +
  geom_point(aes(color = Category), alpha = 0.6, size = 1.2) +
  scale_color_manual(values = cat_colors,
    labels = c("Concordant UP","Concordant DOWN","Discordant","Not significant")) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  geom_text_repel(data = filter(merged_df, symbol %in% all_labels),
    aes(label = symbol), size = 3.2, fontface = "bold",
    box.padding = 0.5, max.overlaps = Inf) +
  labs(title    = "B04 - RNA-seq vs ATAC-seq concordance",
       subtitle = sprintf("R = %.3f | Concordant: %d genes",
                          cor_res$estimate, nrow(concordant)),
       x = "ATAC-seq Log2 Fold Change",
       y = "RNA-seq Log2 Fold Change", color = "") +
  theme_bw() + theme(aspect.ratio = 1, legend.position = "bottom")

save_plot(p_quad, "B04_integration_quadrant", width = 8, height = 8)

# Ranked concordant genes
p_ranked <- concordant %>% head(30) %>%
  mutate(symbol = factor(symbol, levels = rev(symbol))) %>%
  ggplot(aes(x = logFC_RNA, y = symbol,
             color = Category, size = abs(logFC_ATAC))) +
  geom_point(alpha = 0.85) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  scale_color_manual(values = c(Concordant_UP = "firebrick", Concordant_DOWN = "navy")) +
  scale_size_continuous(name = "|ATAC LFC|", range = c(2, 7)) +
  labs(title = "B04 - Top concordant genes (ranked by RNA LFC)",
       subtitle = "Point size = |ATAC Log2FC|",
       x = "RNA-seq Log2 Fold Change", y = NULL, color = NULL) +
  theme_bw()

save_plot(p_ranked, "B04_concordant_ranked", width = 7, height = 8)

# Wnt/PCP pathway highlight
wnt_df <- merged_df %>% filter(symbol %in% wnt_pcp) %>%
  mutate(sig_label = case_when(
    padj_RNA < P_CUTOFF & padj_ATAC < P_CUTOFF ~ "Both significant",
    padj_RNA < P_CUTOFF                          ~ "RNA only",
    padj_ATAC < P_CUTOFF                         ~ "ATAC only",
    TRUE                                         ~ "Not significant"))

if (nrow(wnt_df) > 0) {
  p_wnt <- ggplot(wnt_df, aes(logFC_ATAC, logFC_RNA, color = sig_label)) +
    geom_point(size = 4) +
    geom_text_repel(aes(label = symbol), size = 4, fontface = "bold", max.overlaps = Inf) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    scale_color_manual(values = c("Both significant" = "firebrick", "RNA only" = "#e67e22",
                                  "ATAC only" = "#2980b9", "Not significant" = "grey60")) +
    labs(title = "B04 - Wnt/PCP candidate loci",
         x = "ATAC-seq Log2 Fold Change", y = "RNA-seq Log2 Fold Change", color = "") +
    theme_bw() + theme(aspect.ratio = 1, legend.position = "bottom")
  save_plot(p_wnt, "B04_WntPCP_highlight", width = 7, height = 7)
}

# Fisher's exact test for enrichment of DEGs among DA genes
n_degs <- sum(merged_df$sig_RNA, na.rm = TRUE)
n_da_genes <- sum(merged_df$sig_ATAC, na.rm = TRUE)
n_concordant <- nrow(concordant)
n_universe <- nrow(merged_df)

message("DEGs: ", n_degs)
message("DA genes: ", n_da_genes)
message("Concordant: ", n_concordant)
message("Universe: ", n_universe)

# Build contingency table
fisher_table <- matrix(c(
  n_concordant,                          # DEG & DA
  n_degs - n_concordant,                 # DEG but not DA
  n_da_genes - n_concordant,             # DA but not DEG
  n_universe - n_degs - n_da_genes + n_concordant  # neither
), nrow = 2, byrow = TRUE,
dimnames = list(c("DA", "Not_DA"), c("DEG", "Not_DEG")))

fisher_res <- fisher.test(fisher_table, alternative = "greater")

save_session("B04")
message("\n✓ B04 complete.\n")
