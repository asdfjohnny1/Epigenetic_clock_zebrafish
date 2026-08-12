#!/bin/bash -e
# =============================================================================
# PIPELINE STAGE: 04 / Regulatory region mapping (3 of 6)
# =============================================================================
# SLURM array job: labels every CpG as overlapping a gene body (from
# 01_build_regulatory_annotations.sh) or not, per sample. Same pattern as
# 02_map_cpgs_to_promoters.sh - see that file for the submission command.
# =============================================================================
#SBATCH --job-name=Gene_bodies_methylation
#SBATCH --mail-type=ALL
#SBATCH --mail-user=your-email@institution.ac.uk
#SBATCH -p compute-64-512
#SBATCH --time=7-00:00
#SBATCH --mem=100G
#SBATCH --cpus-per-task=2
#SBATCH -o %x-%A-%a.out
#SBATCH -e %x-%A-%a.err

PROJECT_ROOT="${EPICLOCK_ROOT:-$HOME/epigenetic_clock}"
REG_DIR="$PROJECT_ROOT/multitissue/regulatory_regions"

module load bedtools

mapfile -t FILE_LIST < "$REG_DIR/Transposable_Elements/scripts/data_file_cov.txt"
NUM_FILES=${#FILE_LIST[@]}

if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    echo "SLURM_ARRAY_TASK_ID is not set!"
    exit 1
fi

i=$((SLURM_ARRAY_TASK_ID))
if [[ $i -ge $NUM_FILES ]]; then
    echo "Error: Task ID $i out of range max: $((NUM_FILES-1))"
    exit 1
fi

INPUT_FILE=${FILE_LIST[$i]}
BASE_NAME=$(basename "$INPUT_FILE")
BASE_NAME="${BASE_NAME%%_*}"

COV_FILE="$INPUT_FILE"
GENEBODY_BED="$REG_DIR/gene_bodies_annotated.bed"
OUTPUT="$REG_DIR/Gene_bodies/results"
mkdir -p "$OUTPUT"

COV_GENEB_FILE="$OUTPUT/${BASE_NAME}_CpGs_in_gene_bodies_cov.txt"
COV_NONGENEB_FILE="$OUTPUT/${BASE_NAME}_CpGs_not_in_gene_bodies.txt"
GENEB_LABELLED="$OUTPUT/${BASE_NAME}_CpGs_gene_bodies_labeled.txt"
NONGENEB_LABELLED="$OUTPUT/${BASE_NAME}_CpGs_nongene_bodies_labeled.txt"
COMBINED="$OUTPUT/${BASE_NAME}_CpGs_gene_bodies_and_nongene_bodies_combined.txt"

# Step 1: CpGs overlapping gene bodies
zcat "$COV_FILE" | bedtools intersect -a - -b "$GENEBODY_BED" -wa -wb > "$COV_GENEB_FILE"

# Step 2: CpGs NOT overlapping gene bodies
zcat "$COV_FILE" | bedtools intersect -a - -b "$GENEBODY_BED" -v > "$COV_NONGENEB_FILE"

# Step 3: Label gene-body CpGs
awk 'BEGIN{OFS="\t"} {print $0, "gene_bodies"}' "$COV_GENEB_FILE" > "$GENEB_LABELLED"

# Step 4: Label non-gene-body CpGs, pad with 'NA's for missing gene-body fields
awk 'BEGIN{OFS="\t"} {print $0, "NA", "NA", "NA", "NA", "NA", "NA", "NA", "gene_bodies_no"}' "$COV_NONGENEB_FILE" > "$NONGENEB_LABELLED"

# Step 5: Combine labeled files
cat "$GENEB_LABELLED" "$NONGENEB_LABELLED" > "$COMBINED"

# Step 6: Add sample ID column
awk -v sid="$BASE_NAME" 'BEGIN{OFS="\t"} {print $0, sid}' "$COMBINED" > "${COMBINED%.txt}_tagged.txt"
