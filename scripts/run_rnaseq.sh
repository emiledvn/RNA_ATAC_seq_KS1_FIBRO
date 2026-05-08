#!/usr/bin/env bash
# nf-core/rnaseq 3.14.0 — Kabuki Fibroblast Project
# Annotation: GENCODE v46 / GRCh38
#
# Usage: bash run_rnaseq.sh
# Requirements: 
#   - Nextflow 
#   - Apptainer/Singularity 

set -euo pipefail

PROJECT="${HOME}/analysis/KABUKI_FIBRO_PROJECT"
ANALYSIS_DIR="${PROJECT}/pipeline/RNA"
SAMPLESHEET="${PROJECT}/metadata/samplesheet_RNA.csv"
OUTDIR="${ANALYSIS_DIR}/results"
WORKDIR="/storage/volume01/emile/work/kabuki_rnaseq"
LOGDIR="${ANALYSIS_DIR}/logs"
NXF_CONFIG="${PROJECT}/nextflow.config"
SCREEN_SESSION="kabuki_rnaseq"



# Create output directories
mkdir -p "${OUTDIR}" "${WORKDIR}" "${LOGDIR}"

echo " Kabuki Fibroblast — RNA-seq"
echo " nf-core/rnaseq 3.14.0 | GENCODE v46 | GRCh38"

# Run nf-core/rnaseq pipeline
nextflow run nf-core/rnaseq \
    -revision 3.14.0 \
    -profile apptainer \
    -c "${NXF_CONFIG}" \
    -w "${WORKDIR}" \
    -with-tower \
    -resume \
    --input        "${SAMPLESHEET}" \
    --outdir       "${OUTDIR}" \
    --fasta        'https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/GRCh38.primary_assembly.genome.fa.gz' \
    --gtf          'https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/gencode.v46.annotation.gtf.gz' \
    --save_reference \
    --save_align_intermeds \
    2>&1 | tee "${LOGDIR}/rnaseq_$(date +%Y%m%d_%H%M%S).log"

echo "Done: $(date)"
