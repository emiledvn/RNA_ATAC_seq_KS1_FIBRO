#!/usr/bin/env bash
# nf-core/atacseq 2.1.2 — Kabuki Fibroblast Project
# Reference: GRCh38 via iGenomes
# Broad peak calling (--macs_gsize 2.7e9) suited for chromatin
#
# Usage: bash run_atacseq.sh
# Requirements:
#   - Nextflow 
#   - Apptainer/Singularity 
#   - Conda environment 'nf-run' with nextflow

set -euo pipefail

PROJECT="${HOME}/analysis/KABUKI_FIBRO_PROJECT"
ANALYSIS_DIR="${PROJECT}/pipeline/ATAC"
SAMPLESHEET="${PROJECT}/metadata/samplesheet_ATAC.csv"
OUTDIR="${ANALYSIS_DIR}/results"
WORKDIR="/storage/volume01/emile/work/kabuki_atacseq"
LOGDIR="${ANALYSIS_DIR}/logs"
NXF_CONFIG="${PROJECT}/nextflow.config"
SCREEN_SESSION="kabuki_atacseq"




# Create output directories
mkdir -p "${OUTDIR}" "${WORKDIR}" "${LOGDIR}"

echo " Kabuki Fibroblast — ATAC-seq"
echo " nf-core/atacseq 2.1.2 | GRCh38 | broad peaks"

# Run nf-core/atacseq pipeline
nextflow run nf-core/atacseq \
    -revision 2.1.2 \
    -profile apptainer \
    -c "${NXF_CONFIG}" \
    -w "${WORKDIR}" \
    -with-tower \
    -resume \
    --input          "${SAMPLESHEET}" \
    --outdir         "${OUTDIR}" \
    --genome         GRCh38 \
    --igenomes_base  's3://ngi-igenomes/igenomes' \
    --macs_gsize     2700000000 \
    --save_reference \
    --save_unaligned \
    2>&1 | tee "${LOGDIR}/atacseq_$(date +%Y%m%d_%H%M%S).log"

echo "Done: $(date)"
