#!/bin/bash
# =============================================================================
# PIPELINE STAGE: 04 / Regulatory region mapping (5 of 6)
# =============================================================================
# Run once all three SLURM array jobs (02, 03, 04) have finished. Concatenates
# every per-sample tagged output into one combined file per region type, then
# removes the per-sample intermediates. Run on an interactive/login node, not
# as a batch job.
# =============================================================================
set -euo pipefail
# extglob must be enabled as its own statement, before bash parses any line
# using the !(pattern) syntax below - chaining it with && on the same line
# as the rm command it enables does NOT work reliably.
shopt -s extglob
PROJECT_ROOT="${EPICLOCK_ROOT:-$HOME/epigenetic_clock}"
REG_DIR="$PROJECT_ROOT/multitissue/regulatory_regions"

# --- Promoters ---------------------------------------------------------------
OUTPUT_DIR="$REG_DIR/Promoters/results"
FINAL_PROMO="$OUTPUT_DIR/CpGs_promoters_and_nonpromoters_combined_all_samples.txt"
echo -e "chr\tstart\tend\tmeth_pct\tmeth\tunmeth\tpromoter_chr\tpromoter_start\tpromoter_end\tpromoter_name\tscore\tstrand\tlabel\t(some_promoter_extra_cols)\tsanger_new_id" > "$FINAL_PROMO"
find "$OUTPUT_DIR" -name "*_CpGs_promoters_and_nonpromoters_combined_tagged.txt" -exec cat {} + >> "$FINAL_PROMO"
cd "$OUTPUT_DIR" && rm -vf !("CpGs_promoters_and_nonpromoters_combined_all_samples.txt")

# --- Gene bodies ---------------------------------------------------------------
OUTPUT_DIR="$REG_DIR/Gene_bodies/results"
FINAL_GENEBODY="$OUTPUT_DIR/CpGs_genebodies_and_non_combined_all_samples.txt"
echo -e "chr\tstart\tend\tmeth_pct\tmeth\tunmeth\tgenebody_chr\tgenebody_start\tgenebody_end\tgenebody_name\tscore\tstrand\tlabel\t(some_genebody_extra_cols)\tsanger_new_id" > "$FINAL_GENEBODY"
find "$OUTPUT_DIR" -name "*_CpGs_gene_bodies_and_nongene_bodies_combined_tagged.txt" -exec cat {} + >> "$FINAL_GENEBODY"
cd "$OUTPUT_DIR" && rm -vf !("CpGs_genebodies_and_non_combined_all_samples.txt")

# --- Transposable elements ----------------------------------------------------
OUTPUT_DIR="$REG_DIR/Transposable_Elements/Results"
FINAL_TE="$OUTPUT_DIR/CpGs_TE_and_nonTE_combined_all_samples.txt"
echo -e "chr\tstart\tend\tmeth_pct\tmeth\tunmeth\tte_chr\tte_start\tte_end\tte_name\tscore\tstrand\tsubfamily\tclass\tlabel\tsanger_new_id" > "$FINAL_TE"
find "$OUTPUT_DIR" -name "*_CpGs_TE_and_nonTE_combined_tagged.txt" -exec cat {} + >> "$FINAL_TE"
cd "$OUTPUT_DIR" && rm -vf !("CpGs_TE_and_nonTE_combined_all_samples.txt")

echo "Done. Combined files:"
echo "  $REG_DIR/Promoters/results/CpGs_promoters_and_nonpromoters_combined_all_samples.txt"
echo "  $REG_DIR/Gene_bodies/results/CpGs_genebodies_and_non_combined_all_samples.txt"
echo "  $REG_DIR/Transposable_Elements/Results/CpGs_TE_and_nonTE_combined_all_samples.txt"
echo "These three files feed 06_regulatory_region_models.Rmd."
