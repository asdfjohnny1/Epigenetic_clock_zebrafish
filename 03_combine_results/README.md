# 03 - Combine results

Collates the ~4,320 per-model `.rds` outputs from stage 02 into the master
summary tables used by everything downstream. Run once - this is the single
shared results-aggregation step for **both aims** of the study.

| Script | What it does |
|---|---|
| `01_collate_model_metrics.R` | Reads every per-model `.rds` file; extracts predicted-vs-chronological ages (train/test/validation), MAE, R^2, and selected-CpG sets; incrementally writes three collated CSVs (metrics, predicted ages, CpG sets). |
| `02_build_summary_tables.R` | Reads the CSVs from the previous script (plus the raw `.rds` files again, for CpG identities); builds the final summary tables and diagnostic scatterplots (predicted vs. chronological age per tissue/model). |

Run in order: `01` then `02`.

**3,036 of the 4,320 trained models pass QC/performance criteria and are
retained** for all analyses reported in the manuscript (see Methods).

Output consumed by both aim-specific analysis notebook sets:
[`05_clock_design_benchmarking/`](../05_clock_design_benchmarking/) and
[`06_cross_tissue_ageing_biology/`](../06_cross_tissue_ageing_biology/).
