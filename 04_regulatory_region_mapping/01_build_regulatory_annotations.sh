#!/bin/bash
# =============================================================================
# PIPELINE STAGE: 04 / Regulatory region mapping (1 of 6)
# =============================================================================
# Downloads the Ensembl GTF annotation for the zebrafish genome (GRCz11) and
# derives two BED files used by the mapping scripts that follow:
#   - promoters_annotated.bed   (+/-1 kb around each transcript's TSS)
#   - gene_bodies_annotated.bed (full gene span, from the "gene" GTF features)
# Run once, before 02_map_cpgs_to_promoters.sh / 03_map_cpgs_to_gene_bodies.sh.
# =============================================================================
set -euo pipefail

# Set the root directory of this project (edit here, or export EPICLOCK_ROOT
# before running). All paths below are built from this.
PROJECT_ROOT="${EPICLOCK_ROOT:-$HOME/epigenetic_clock}"
REG_DIR="$PROJECT_ROOT/multitissue/regulatory_regions"
mkdir -p "$REG_DIR"
cd "$REG_DIR"

# Download regulatory features GTF file
wget ftp://ftp.ensembl.org/pub/release-111/gtf/danio_rerio/Danio_rerio.GRCz11.111.gtf.gz
gunzip Danio_rerio.GRCz11.111.gtf.gz

# --- Promoter regions -------------------------------------------------------
# Each row corresponds to a putative promoter region (TSS +/- 1 kb) for a
# specific transcript, derived from the GTF annotation.
gawk '$3 == "transcript"' Danio_rerio.GRCz11.111.gtf | \
gawk 'BEGIN{OFS="\t"} {
    match($0, /gene_id "[^"]+"/, gid);
    match($0, /transcript_id "[^"]+"/, tid);
    match($0, /gene_name "[^"]+"/, gname);
    gene_id = (gid[0] != "") ? gid[0] : "gene_id \"NA\"";
    transcript_id = (tid[0] != "") ? tid[0] : "transcript_id \"NA\"";
    gene_name = (gname[0] != "") ? gname[0] : "gene_name \"NA\"";

    if ($7 == "+") {
        start = $4 - 1000; end = $4 + 1000;
    } else {
        start = $5 - 1000; end = $5 + 1000;
    }
    if (start < 0) start = 0;
    print $1, start, end, gene_name, ".", $7, gene_id, transcript_id;
}' > promoters_annotated.bed

# --- Gene body regions -------------------------------------------------------
# Full transcribed span of each gene (promoter through 3' end).
gawk '$3 == "gene"' Danio_rerio.GRCz11.111.gtf | \
gawk 'BEGIN{OFS="\t"}
{
    match($0, /gene_id "[^"]+"/, gid);
    match($0, /gene_name "[^"]+"/, gname);
    gene_id = (gid[0] != "") ? gid[0] : "gene_id \"NA\"";
    gene_name = (gname[0] != "") ? gname[0] : "gene_name \"NA\"";
    start = $4 - 1;  # BED is 0-based
    end = $5;        # BED is half-open [start, end)
    print $1, start, end, gene_name, ".", $7, gene_id;
}' > gene_bodies_annotated.bed

echo "Wrote $REG_DIR/promoters_annotated.bed and $REG_DIR/gene_bodies_annotated.bed"

# -----------------------------------------------------------------------------
# Why promoters and gene bodies?
# Promoters are regulatory DNA sequences just upstream of a gene's
# transcription start site, controlling initiation of transcription.
# The gene body spans from just after the promoter/TSS through all exons and
# introns to the transcription termination site; methylation and mutation
# patterns there can influence gene expression, RNA processing and genomic
# stability. Promoter methylation is often associated with gene silencing;
# gene body methylation is linked to active transcriptional regulation.
# Differences in these patterns can indicate regions of regulatory change,
# biological ageing, or environmental/stress responses.
# -----------------------------------------------------------------------------
