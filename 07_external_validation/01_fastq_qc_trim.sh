#!/bin/bash
# ==========================================================================
# PIPELINE STAGE: 07 / External validation - step 1: FastQ QC + trimming
# ==========================================================================
# Runs Trim Galore on one single-end FastQ sample from the external
# validation dataset (Mayne et al. 2020's public zebrafish RRBS dataset -
# see 00_download_external_dataset.sh - is single-end). Called once per
# SLURM array task by hpc_array_jobs/submit_01_fastq_qc_trim.sh, which
# supplies the sample ID as $1.
#
# Uses --rrbs (removes the 2 bp MspI fill-in artefact from RRBS library
# prep), matching how both this study's own data and the external dataset
# above were generated. Drop that flag if you point this at a WGBS
# dataset instead.
#
# If you point this pipeline at a different, PAIRED-END external dataset,
# replace the single trim_galore call below with:
#   trim_galore --rrbs --paired --cores 4 --output_dir "$OUT_DIR" \
#     "$RAW_DIR/${SAMPLE}_R1.fastq.gz" "$RAW_DIR/${SAMPLE}_R2.fastq.gz"
# (output becomes <sample>_R1_val_1.fq.gz / <sample>_R2_val_2.fq.gz - also
# update the -1/-2 arguments in 02_bismark_align_extract.sh accordingly).
#
# Input:  $PROJECT_ROOT/external_validation/FastQ_raw/<sample>.fastq.gz
# Output: $PROJECT_ROOT/external_validation/FastQ_trimmed/<sample>_trimmed.fq.gz
# ==========================================================================

set -euo pipefail

PROJECT_ROOT="${EPICLOCK_ROOT:-$(pwd)}"
SAMPLE="$1"

RAW_DIR="$PROJECT_ROOT/external_validation/FastQ_raw"
OUT_DIR="$PROJECT_ROOT/external_validation/FastQ_trimmed"
mkdir -p "$OUT_DIR"

echo "[$(date)] Trimming sample: $SAMPLE"

trim_galore \
  --rrbs \
  --cores 4 \
  --output_dir "$OUT_DIR" \
  "$RAW_DIR/${SAMPLE}.fastq.gz"

echo "[$(date)] Done: $SAMPLE"
