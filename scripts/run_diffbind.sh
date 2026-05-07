#!/usr/bin/env bash
# DiffBind analysis with screen session management

set -euo pipefail

PROJECT="${HOME}/analysis/KABUKI_FIBRO_PROJECT"
SCRIPT_DIR="${PROJECT}/scripts"
LOGDIR="${PROJECT}/results/ATAC/DiffBind/logs"
SCREEN_SESSION="diffbind"

# Launch in screen if not already running in one
if [[ -z "${STY:-}" ]]; then
    mkdir -p "${LOGDIR}"
    exec screen -S "${SCREEN_SESSION}" bash "$0" "$@"
fi

# Activate conda environment
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate atac_r

mkdir -p "${LOGDIR}"

# Run DiffBind analysis with logging
Rscript "${SCRIPT_DIR}/diffbind_analysis.R" \
    2>&1 | tee "${LOGDIR}/diffbind_$(date +%Y%m%d_%H%M%S).log"
