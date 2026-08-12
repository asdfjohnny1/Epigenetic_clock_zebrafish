#!/bin/bash -e
# =============================================================================
# PIPELINE STAGE: 04 / Regulatory region mapping (2 of 6)
# =============================================================================
# SLURM array job: for each sample's methylation coverage file (.cov.gz),
# labels every CpG as overlapping a promoter (from 01_build_regulatory_annotations.sh)
# or not, using bedtools, and writes one combined+tagged output file per sample.
#
# Submit with:
#   NUM_FILES=$(wc -l < $PROJECT_ROOT/multitissue/regulatory_regions/Transposable_Elements/scripts/data_file_cov.txt)
#   sbatch --array=0-$((NUM_FILES - 1))%100 02_map_cpgs_to_promoters.sh
# (the same sample file list, produced in 04_map_cpgs_to_transposable_elements.sh,
# is reused here so all three region types are mapped over the same samples)
# =============================================================================
#SBATCH --job-name=Promoter_methylation
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
PROMOTER_BED="$REG_DIR/promoters_annotated.bed"
OUTPUT="$REG_DIR/Promoters/results"
mkdir -p "$OUTPUT"

COV_PROMO_FILE="$OUTPUT/${BASE_NAME}_CpGs_in_promoters_cov.txt"
COV_NONPROMO_FILE="$OUTPUT/${BASE_NAME}_CpGs_not_in_promoters.txt"
PROMO_LABELLED="$OUTPUT/${BASE_NAME}_CpGs_promoters_labeled.txt"
NONPROMO_LABELLED="$OUTPUT/${BASE_NAME}_CpGs_nonpromoters_labeled.txt"
COMBINED="$OUTPUT/${BASE_NAME}_CpGs_promoters_and_nonpromoters_combined.txt"

# Step 1: CpGs overlapping promoters
zcat "$COV_FILE" | bedtools intersect -a - -b "$PROMOTER_BED" -wa -wb > "$COV_PROMO_FILE"

# Step 2: CpGs NOT overlapping promoters
zcat "$COV_FILE" | bedtools intersect -a - -b "$PROMOTER_BED" -v > "$COV_NONPROMO_FILE"

# Step 3: Label promoter CpGs
awk 'BEGIN{OFS="\t"} {print $0, "promoter"}' "$COV_PROMO_FILE" > "$PROMO_LABELLED"

# Step 4: Label non-promoter CpGs, pad with 'NA's for missing promoter fields
awk 'BEGIN{OFS="\t"} {print $0, "NA", "NA", "NA", "NA", "NA", "NA","NA", "NA","promoter_no"}' "$COV_NONPROMO_FILE" > "$NONPROMO_LABELLED"

# Step 5: Combine labeled files
cat "$PROMO_LABELLED" "$NONPROMO_LABELLED" > "$COMBINED"

# Step 6: Add sample ID column
awk -v sid="$BASE_NAME" 'BEGIN{OFS="\t"} {print $0, sid}' "$COMBINED" > "${COMBINED%.txt}_tagged.txt"
