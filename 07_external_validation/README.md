# 07 - External validation

Applies each tissue's trained clock (the `final_clock` stored in
`data/Metrics_all_models/best_models/best_model_<tissue>.rds`) to an
independent zebrafish methylation dataset, from raw FastQ through to
predicted age and, where chronological age is known for the external
samples, a direct comparison against this study's own validation
performance. Addresses Action item 7.

Defaults to the dataset behind Mayne et al. (2020), "A DNA methylation age
predictor for zebrafish," *Aging* 12:24817-24835
([doi:10.18632/aging.202400](https://doi.org/10.18632/aging.202400)) - 96
single-end RRBS samples of caudal fin tissue, publicly available with no
access request needed. `sample_list.txt` and
`external_metadata_template.csv` are already filled in with its 96 real
sample IDs, tissue, age (in weeks, from the paper's Supplementary Table 1),
and sex.

**Units:** this study's own Age metadata - and therefore every trained
clock's predictions - is in **months** (confirmed range 3-22 in
`data/final_metadata_uniqueID_new.xlsx`'s `multitissue` cohort). The age
column in `external_metadata_template.csv` must be named `age_months`,
`age_weeks`, or `age_days` (unit read from the column name itself, not a
separate field, so it can't be set inconsistently) - `03_apply_best_clocks.R`
converts it to months before comparing against predicted age. Getting this
wrong doesn't corrupt R² (Pearson correlation is scale-invariant) but does
inflate/distort MAE and makes the predicted-vs-actual plot's 1:1 reference
line meaningless.

Run in order:

| # | Script | What it does |
|---|---|---|
| 0 | `00_download_external_dataset.sh` | Downloads all 96 raw FastQ files from the CSIRO Data Access Portal (verified working, no auth needed - see script header for the exact source and its caveats). |
| 1 | `01_fastq_qc_trim.sh` | Trim Galore QC/trimming of one single-end RRBS sample (called per SLURM array task by `hpc_array_jobs/submit_01_fastq_qc_trim.sh`). |
| 2 | `02_bismark_align_extract.sh` | Bismark alignment to danRer11 and methylation extraction for one sample, producing a standard `.bismark.cov.gz` coverage file (called per SLURM array task by `hpc_array_jobs/submit_02_bismark_align_extract.sh`). |
| 3 | `03_apply_best_clocks.R` | Extracts methylation at exactly the CpGs/windows each tissue's clock was trained on, median-imputes any uncovered ones using this study's own training-set medians, predicts age via the stored `glmnet` fit, and reports validation MAE/R² against this study's own numbers wherever chronological age is known (single job, `hpc_array_jobs/submit_03_apply_best_clocks.sh`). |
| 4 | `04_diagnose_methylation_shift.R` | Only needed if step 3's validation looks poor. Re-extracts the external samples' raw methylation at the clock's CpGs and compares per-CpG mean/SD against the training data at the same windows (batch/platform-shift check), then separately computes each CpG's own correlation with age in the training cohort vs. the external cohort (checks whether these specific CpGs carry any age signal at all in an independent cohort, as opposed to just being shifted) (single job, `hpc_array_jobs/submit_04_diagnose_methylation_shift.sh`). |

## Setup

1. Download the danRer11 genome FASTA (UCSC or Ensembl) into
   `$EPICLOCK_ROOT/external_validation/genome/` and run
   `bismark_genome_preparation` on it once.
2. Run `bash 00_download_external_dataset.sh` (needs `curl` and `python3`;
   ~130 GB total across 96 files, so make sure there's disk space and
   budget real time for it). It writes `sample_list.txt` itself from the
   real file listing.
3. Chronological age (`age_weeks`) and sex are already filled in for all
   96 samples in `external_metadata_template.csv`, sourced by hand from
   the paper's Supplementary Table 1 (`aging-12-202400-s002.xlsx` - both
   the journal site and PMC block scripted downloads of it, so this had
   to be done manually via
   [PMC7803548](https://pmc.ncbi.nlm.nih.gov/articles/PMC7803548/)). If
   you point this pipeline at a different dataset, name its age column
   `age_months`, `age_weeks`, or `age_days` to match.

```bash
export EPICLOCK_ROOT=/path/to/your/project/root
bash 07_external_validation/00_download_external_dataset.sh
sbatch hpc_array_jobs/submit_01_fastq_qc_trim.sh
sbatch hpc_array_jobs/submit_02_bismark_align_extract.sh
sbatch hpc_array_jobs/submit_03_apply_best_clocks.sh
```

Output: `results/external_validation/predicted_ages.csv`, plus
`validation_summary.csv` and `predicted_vs_actual.png` for any tissue with
at least one sample of known age.

**First run result:** on this dataset, the fin clock does not transfer -
predicted vs. actual age R² = 0.005 (vs. this study's own 0.83), with
~11% of samples predicted a negative age (clamped to 0). Ruled out as
causes: CpG coverage is good (median 20/22), there's no chromosome-naming
mismatch (`GRCz11`/Ensembl and `danRer11`/UCSC are the same assembly and
coordinates, just `chr`-prefixed or not), and `04_diagnose_methylation_
shift.R` shows no broad batch/platform shift either - 14 of 22 CpGs match
the training data's mean methylation within 3 percentage points (2 clear
outliers aside: `chr20:51030282` and `chr21:9782310`).

Most likely explanation instead: these 22 CpGs were selected by elastic
net because they correlate with age *in this study's own cohort* - there
is no reason to expect the same loci to carry any age signal in an
independent cohort collected for a *different* clock study (Mayne et
al.'s own model uses a different, non-overlapping 26-29 CpGs). The same
diagnostic script's second half tests this directly: per-CpG correlation
with age, training cohort vs. external cohort (`age_correlation_
diagnostic.png` / the `train_age_cor` / `external_age_cor` columns in
`methylation_shift_diagnostic.csv`). If CpGs that are clearly
age-correlated in training show ~zero correlation with age externally,
that confirms it - and is itself a legitimate, reportable limitation for
Action item 7 (narrow, elastic-net-selected CpG panels are not expected
to be portable across independently-ascertained cohorts) rather than
something to "fix."

`best_model_all_organs.rds` (the `Multi_organ` clock) is excluded from
this repository - see [`data/README.md`](../data/README.md) - so
validating against it requires regenerating that file via
[`03_combine_results/`](../03_combine_results/) first. This particular
external dataset is caudal fin only, so it can only validate
`best_model_fin.rds` regardless.

To point this pipeline at a *different* external dataset instead: replace
step 0 with your own download step, and see the "PAIRED-END" comment
blocks in `01_fastq_qc_trim.sh` and `02_bismark_align_extract.sh` if it's
paired-end rather than single-end (both scripts currently default to
single-end, matching the Mayne et al. dataset above). Both also assume
RRBS library prep (`--rrbs` trimming, no deduplication), matching this
study's own data; adjust if the new dataset is WGBS instead (see each
script's header comment).
