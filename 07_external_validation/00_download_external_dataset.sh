#!/bin/bash
# ==========================================================================
# PIPELINE STAGE: 07 / External validation - step 0: download raw data
# ==========================================================================
# Downloads the raw RRBS FastQ data behind Mayne et al. (2020), "A DNA
# methylation age predictor for zebrafish," Aging 12:24817-24835,
# https://doi.org/10.18632/aging.202400 - 96 single-end RRBS samples from
# caudal fin tissue (AB strain), the dataset flagged in Action item 7 as
# the obvious public zebrafish DNAm dataset to validate against.
#
# Data source (verified working as of this writing): CSIRO Data Access
# Portal collection "Reduced representation bisulfite sequencing of
# caudal fin Zebrafish", https://doi.org/10.25919/5f63ce026960a
# (DAP collection ID 50185). Its file-listing API returns time-limited
# presigned S3 URLs (~48h expiry), which is why this script re-fetches
# the listing itself each run rather than hardcoding URLs.
#
# NOTE: this is single-end sequencing (one FastQ file per sample, no R2) -
# different from this study's own paired-end data. 01_fastq_qc_trim.sh
# and 02_bismark_align_extract.sh have been written accordingly; if you
# point this pipeline at a different, paired-end external dataset instead,
# see the "paired-end" comment blocks in those two scripts to switch back.
#
# Chronological age per sample is NOT included in this file listing - it
# lives in the paper's Supplementary Table 1 (aging-12-202400-s002.xlsx).
# Both the journal site and PMC block scripted downloads of it (bot
# challenges); download it by hand from
# https://pmc.ncbi.nlm.nih.gov/articles/PMC7803548/ (Supplementary
# Materials section) and use it to fill in the `age` column of
# external_metadata_template.csv after this script runs.
#
# Output: $PROJECT_ROOT/external_validation/FastQ_raw/<sample>.fastq.gz
# ==========================================================================

set -euo pipefail

PROJECT_ROOT="${EPICLOCK_ROOT:-$(pwd)}"
RAW_DIR="$PROJECT_ROOT/external_validation/FastQ_raw"
mkdir -p "$RAW_DIR"

DAP_COLLECTION_ID=50185
DAP_API="https://data.csiro.au/dap/ws/v2/collections/${DAP_COLLECTION_ID}/data"
UA="Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36"

echo "[$(date)] Fetching file listing from CSIRO Data Access Portal..."
LISTING_JSON="$RAW_DIR/.dap_listing.json"
curl -sL -A "$UA" -H "Accept: application/json" "$DAP_API" -o "$LISTING_JSON"

N_FILES=$(python3 -c "import json; print(len(json.load(open('$LISTING_JSON'))['file']))")
echo "Found $N_FILES files in the collection."

# -----------------------
# Download every file, skipping any already present (safe to re-run;
# re-fetch the listing above first if a previous run's presigned URLs
# have expired)
# -----------------------
python3 - "$LISTING_JSON" "$RAW_DIR" <<'PYEOF'
import json, sys, subprocess, pathlib

listing_path, raw_dir = sys.argv[1], sys.argv[2]
with open(listing_path) as f:
    files = json.load(f)["file"]

sample_list = []
for entry in sorted(files, key=lambda f: f["filename"]):
    filename = entry["filename"]
    sample_id = filename.split("_")[0]
    sample_list.append(sample_id)

    dest = pathlib.Path(raw_dir) / f"{sample_id}.fastq.gz"
    if dest.exists() and dest.stat().st_size > 0:
        print(f"  [skip] {sample_id} (already downloaded)")
        continue

    url = entry["presignedLink"]["href"]
    print(f"  [get]  {sample_id}  <-  {filename}  ({entry['fileSize'] / 1e9:.2f} GB)")
    subprocess.run(["curl", "-sL", "--fail", "-o", str(dest), url], check=True)

# Write sample_list.txt alongside this script for convenience - review/
# edit before submitting the array jobs (e.g. to restrict to a subset).
out_list = pathlib.Path(raw_dir).parent.parent / "07_external_validation" / "sample_list.txt"
out_list.write_text("\n".join(sample_list) + "\n")
print(f"\nWrote {len(sample_list)} sample IDs to {out_list}")
PYEOF

rm -f "$LISTING_JSON"

echo ""
echo "[$(date)] Done. FastQ files are in: $RAW_DIR"
echo "Next: fill in the 'age' column of external_metadata_template.csv"
echo "using Supplementary Table 1 from the paper (see header comment"
echo "above for the manual-download link), then continue to"
echo "01_fastq_qc_trim.sh."
