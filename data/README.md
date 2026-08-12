# Data

Most small/derived datasets needed to run the downstream analysis notebooks
(stages 05/06) are included directly in this folder. Large raw/intermediate
files are not - see "Excluded from this repository" below.

## Included

```
data/
├── final_metadata_uniqueID_new.xlsx           sample metadata (age, sex, tissue, batch, Sanger IDs)
├── Bracken_et_al_2006_PCGT_list/
│   ├── Bracken2006_PCCT_human.csv              human Polycomb Group Target (PCGT) gene list
│   └── BrackenSuppTable3.pdf                   source supplementary table (Bracken et al. 2006)
├── Metrics_all_models/
│   ├── final_merged/
│   │   ├── results_metrics.csv / _final.csv    MAE/R^2 per model (output of 03_combine_results/)
│   │   ├── results_predicted.csv / _final.csv  predicted vs. chronological age per model
│   │   └── results_cpg.csv / _final.csv        selected CpGs per model
│   ├── best_models/
│   │   ├── best_model_brain.rds
│   │   ├── best_model_fin.rds
│   │   └── best_model_intestines.rds
│   └── methylation_data/
│       └── methylation_cpgs.csv                methylation values for the shared/overlapping CpG set
└── Regulatory_Regions/
    └── RR.rds                                   summarised regulatory-region mapping (small; used by 04_regulatory_region_mapping/06_regulatory_region_models.Rmd)
```

Note: `_final` suffixed files are a later-refined version of their
non-suffixed counterpart, both kept because different notebook parts read
one or the other (check the specific `.Rmd` if you need to know which).

## Excluded from this repository

These are either too large for a normal git repository (GitHub hard-blocks
any single file over 100MB) or are raw sequencing data under ENA embargo.
Regenerate them by running the earlier pipeline stages, or source them
separately (Git LFS / Zenodo / OSF companion deposit - your call, not set up
here):

| Path (if present) | Size | Why excluded | How to get it |
|---|---|---|---|
| `Metrics_all_models/best_models/best_model_all_organs.rds` | ~254 MB | Exceeds GitHub's 100 MB hard limit | Regenerate via [`03_combine_results/`](../03_combine_results/) |
| `Metrics_all_models/methylation_data/*.tsv` | 0.5-3.8 GB each | Far exceeds 100 MB limit | Regenerate via [`01_annotation/`](../01_annotation/) |
| `Regulatory_Regions/final_RR_df1.csv` | ~250 MB | Exceeds 100 MB limit | Regenerate via [`04_regulatory_region_mapping/`](../04_regulatory_region_mapping/) (needed by `06_cross_tissue_ageing_biology/09_regulatory_regions_analysis.Rmd`) |
| `FastQ_files/` | ~15 GB | Raw Bismark coverage files - this **is** the ENA-deposited raw data | Download from ENA accession **PRJEB107117** once the embargo lifts (see below) |

**Practical consequence:** stages 01, 02 (which operate on raw coverage
files) and one specific part of stage 06
(`09_regulatory_regions_analysis.Rmd`) cannot be run end-to-end from this
repository alone. Everything from stage 03 onward that depends only on the
files listed as "Included" above will run.

**Raw RRBS data:** ENA accession **PRJEB107117**
(https://www.ebi.ac.uk/ena/browser/view/PRJEB107117), under embargo until
acceptance for publication.

## Expected local layout for stages 01-04

The shell/HPC pipeline stages (01, 02, 04) assume a separate `PROJECT_ROOT`
directory (see each stage's README / script header - set via the
`EPICLOCK_ROOT` environment variable), distinct from this `data/` folder,
containing:

```
$PROJECT_ROOT/
├── multitissue/
│   ├── making_annotation_files/     # stage 01 output
│   ├── making_clocks/R/results*/    # stage 02 output (per-tissue, per-CpG-grouping)
│   └── regulatory_regions/          # stage 04 output
└── cov_files/                       # Bismark .cov.gz coverage files (input to stage 01)
```

The aim-specific analysis notebooks (`05_clock_design_benchmarking/`,
`06_cross_tissue_ageing_biology/`, and `04_regulatory_region_mapping/06_regulatory_region_models.Rmd`)
instead read from *this* `data/` folder via the `here` package (e.g.
`here("data/Metrics_all_models/best_models/...")`) - that's why its contents
are tracked directly rather than living under `$PROJECT_ROOT`.
