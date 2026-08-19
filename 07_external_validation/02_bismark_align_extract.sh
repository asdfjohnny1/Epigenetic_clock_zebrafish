#!/bin/bash
# ==========================================================================
# PIPELINE STAGE: 07 / External validation - step 2: Bismark alignment + methylation extraction
# ==========================================================================
# Aligns one trimmed single-end sample to the zebrafish genome (danRer11)
# with Bismark, and extracts per-CpG methylation calls to a standard
# Bismark .cov.gz coverage file - the same format
# 01_annotation/making_annotation.R reads for this study's own data, so
# 03_apply_best_clocks.R can extract methylation at the trained clocks'
# CpGs directly from this output. Called once per SLURM array task by
# hpc_array_jobs/submit_02_bismark_align_extract.sh, which supplies the
# sample ID as $1.
#
# One-time setup before running this script:
#   1. Download the danRer11 genome FASTA (UCSC or Ensembl) into
#      $PROJECT_ROOT/external_validation/genome/
#   2. bismark_genome_preparation $PROJECT_ROOT/external_validation/genome/
#
# Deduplication is skipped by default: like this study's own data, RRBS
# fragments are size-selected rather than randomly sheared, so PCR
# duplicates cannot be identified by alignment position alone. Uncomment
# the deduplicate_bismark block below if you point this at a WGBS dataset.
#
# If you point this pipeline at a different, PAIRED-END external dataset,
# replace -1/-2 below with paired inputs and add --paired-end to the
# methylation-extractor call (see the equivalent note in
# 01_fastq_qc_trim.sh for the matching trim_galore change).
#
# Input:  $PROJECT_ROOT/external_validation/FastQ_trimmed/<sample>_trimmed.fq.gz
# Output: $PROJECT_ROOT/external_validation/cov_files/<sample>.bismark.cov.gz
# ==========================================================================

set -euo pipefail

PROJECT_ROOT="${EPICLOCK_ROOT:-$(pwd)}"
SAMPLE="$1"

TRIM_DIR="$PROJECT_ROOT/external_validation/FastQ_trimmed"
GENOME_DIR="$PROJECT_ROOT/external_validation/genome"
ALIGN_DIR="$PROJECT_ROOT/external_validation/bismark_aligned"
COV_DIR="$PROJECT_ROOT/external_validation/cov_files"
mkdir -p "$ALIGN_DIR" "$COV_DIR"

echo "[$(date)] Aligning sample: $SAMPLE"

bismark \
  --genome "$GENOME_DIR" \
  "$TRIM_DIR/${SAMPLE}_trimmed.fq.gz" \
  --output_dir "$ALIGN_DIR" \
  --multicore 2

BAM="$ALIGN_DIR/${SAMPLE}_trimmed_bismark_bt2.bam"

# --- WGBS only: deduplicate before methylation extraction ---
# deduplicate_bismark --single-end --bam "$BAM" --output_dir "$ALIGN_DIR"
# BAM="$ALIGN_DIR/$(basename "${BAM%.bam}").deduplicated.bam"

echo "[$(date)] Extracting methylation calls: $SAMPLE"

bismark_methylation_extractor \
  --single-end \
  --comprehensive \
  --bedGraph \
  --genome_folder "$GENOME_DIR" \
  --output "$ALIGN_DIR" \
  "$BAM"

# bismark_methylation_extractor writes <sample>...bismark.cov.gz into
# $ALIGN_DIR as part of --bedGraph; copy it into cov_files/ under a plain
# <sample>.bismark.cov.gz name so downstream sample-ID matching is simple.
COV_OUT=$(ls "$ALIGN_DIR/${SAMPLE}"*.bismark.cov.gz | head -n1)
cp "$COV_OUT" "$COV_DIR/${SAMPLE}.bismark.cov.gz"

echo "[$(date)] Done: $SAMPLE"
