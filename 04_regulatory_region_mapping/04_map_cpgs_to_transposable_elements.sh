#!/bin/bash -e
# =============================================================================
# PIPELINE STAGE: 04 / Regulatory region mapping (4 of 6)
# =============================================================================
# Builds the shared sample file list (data_file_cov.txt, reused by
# 02_map_cpgs_to_promoters.sh and 03_map_cpgs_to_gene_bodies.sh), then runs a
# SLURM array job labelling every CpG as overlapping a RepeatMasker
# Transposable Element (TE) or not, per sample.
#
# Requires a TE annotation BED file at:
#   $PROJECT_ROOT/multitissue/regulatory_regions/Transposable_Elements/data/temp_TE2.bed
# (RepeatMasker output, filtered to non-redundant autosomal elements with
# Class/subfamily information - see Methods).
# =============================================================================
#SBATCH --job-name=TE_methylation
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
TE_DIR="$REG_DIR/Transposable_Elements"
mkdir -p "$TE_DIR/scripts" "$TE_DIR/Results" "$TE_DIR/log"

# Build the sample file list (run once, outside the array job, before
# submitting this script - or run this script itself once with
# SLURM_ARRAY_TASK_ID unset to just (re)generate the list and exit early).
if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    cd "$PROJECT_ROOT/cov_files/Cross_tissue_clock/All_organs"
    realpath *cov.gz > "$TE_DIR/scripts/data_file_cov.txt"
    echo "Sample list written to $TE_DIR/scripts/data_file_cov.txt"
    echo "Now submit as an array job:"
    echo '  NUM_FILES=$(wc -l < '"$TE_DIR/scripts/data_file_cov.txt"')'
    echo '  sbatch --array=0-$((NUM_FILES - 1))%100 04_map_cpgs_to_transposable_elements.sh'
    exit 0
fi

module load bedtools

mapfile -t FILE_LIST < "$TE_DIR/scripts/data_file_cov.txt"
NUM_FILES=${#FILE_LIST[@]}

i=$((SLURM_ARRAY_TASK_ID))
if [[ $i -ge $NUM_FILES ]]; then
    echo "Error: Task ID $i is out of range (max: $((NUM_FILES - 1)))"
    exit 1
fi

INPUT_FILE=${FILE_LIST[$i]}
BASE_NAME=$(basename "$INPUT_FILE")
BASE_NAME="${BASE_NAME%%_*}"

COV_FILE="$INPUT_FILE"
TE_BED="$TE_DIR/data/temp_TE2.bed"
OUTPUT="$TE_DIR/Results"

COV_TE_FILE="$OUTPUT/${BASE_NAME}_CpGs_in_TEs_cov.txt"
COV_NONTE_FILE="$OUTPUT/${BASE_NAME}_CpGs_not_in_TEs.txt"
TE_LABELLED="$OUTPUT/${BASE_NAME}_CpGs_TE_labeled.txt"
NONTE_LABELLED="$OUTPUT/${BASE_NAME}_CpGs_nonTE_labeled.txt"
COMBINED="$OUTPUT/${BASE_NAME}_CpGs_TE_and_nonTE_combined.txt"

# Step 1: Intersect CpGs with TEs
zcat "$COV_FILE" | bedtools intersect -a - -b "$TE_BED" -wa -wb > "$COV_TE_FILE"

# Step 2: Extract CpGs not overlapping TEs
zcat "$COV_FILE" | bedtools intersect -a - -b "$TE_BED" -v > "$COV_NONTE_FILE"

# Step 3: Add TE label to TE CpGs (13 columns total)
awk 'BEGIN{OFS="\t"} {print $0, "TE"}' "$COV_TE_FILE" > "$TE_LABELLED"

# Step 4: Add dummy TE fields + nonTE label to non-TE CpGs (pad to 15 columns)
awk 'BEGIN{OFS="\t"} {print $0, "NA", "NA", "NA", "NA", "NA", "NA","NA","NA", "nonTE"}' "$COV_NONTE_FILE" > "$NONTE_LABELLED"

# Step 5: Merge labelled TE and nonTE CpGs
cat "$NONTE_LABELLED" "$TE_LABELLED" > "$COMBINED"

# Step 6: Add sample ID as column 'sanger_new_id' to combined file
awk -v sid="$BASE_NAME" 'BEGIN{OFS="\t"} {print $0, sid}' "$COMBINED" > "${COMBINED%.txt}_tagged.txt"
