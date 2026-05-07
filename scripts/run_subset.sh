#!/bin/bash
# Extract specific gene regions from BAM files for visualization

set -e

PROJECT_DIR="$HOME/analysis/ED25_005_all_fibro_RNA"

cd "$PROJECT_DIR"

# Extract gene coordinates from GENCODE GTF
awk -F '\t' -v OFS='\t' '$3 == "gene" && $9 ~ /gene_name "(KMT2D|KDM6A|XIST|KDM5D|RPS4Y1|GAPDH)"/ {print $1, $4-1, $5}' \
  results_gencode/genome/*.gtf > targets.bed

# Append common SNPs for sample fingerprinting (GRCh38)
cat >> targets.bed << SNPS
chr17	7673801	7673803	# rs1042522 (TP53)
chr1	155161622	155161624	# rs321198 (MUC1)
chr19	41354311	41354313	# rs1800470 (TGFB1)
SNPS

mkdir -p subsetted_bams

# Extract regions from each BAM file
for bam in results_gencode/star_salmon/*.markdup.sorted.bam; do
    [ -e "$bam" ] || continue
    filename=$(basename "$bam")
    samtools view -b -L targets.bed "$bam" > "subsetted_bams/$filename"
    samtools index "subsetted_bams/$filename"
done
