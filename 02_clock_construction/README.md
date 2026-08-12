# 02 - Clock construction

Builds all 4,320 candidate epigenetic clock configurations (3 CpG-aggregation
strategies x 2 regression families x 4 feature filters x 5 R^2 thresholds x
up to 12 Tile window sizes x 4 tissue groupings). Two sub-stages:

## `step1-5_build_candidate_clocks/`

Data preparation and probe filtering, per CpG-aggregation strategy:

| Script | CpG aggregation |
|---|---|
| `point_clock.R` | Point (individual CpGs) |
| `tile_clock.R` | Tile (fixed genomic windows, whole-RRBS or CpG-island-restricted) |
| `cluster_clock.R` | Cluster (DBSCAN-based correlated-CpG groups) |

Each script runs steps 1-5: data prep -> probe filtering -> per-probe R^2 /
mean methylation -> modelling-scaffold initialisation -> R^2-threshold
filtering. Outputs `*_clock_setup.rds` scaffold files, one per
tissue x R^2-threshold combination, consumed by step 6.

## `step6_final_model_fitting/`

Step 6: the nested iterative model-fitting procedure described in the
Methods (100 outer x 50 inner iterations; random candidate-CpG subsetting;
held-out validation; best-outer-iteration selection). One R script per
CpG-grouping x regression-family combination (6 total: Point/Tile/Cluster x
Lasso/Elastic Net); each is **array-ready** - one invocation processes one
scaffold file (one SLURM array task).

Submit via the matching script in `hpc_array_jobs/`, e.g.:

```bash
export EPICLOCK_ROOT=/path/to/your/project/root
NUM_FILES=$(wc -l < "$EPICLOCK_ROOT/multitissue/making_clocks/R/results_split_organ/point_clock/step4/scaffold_file_list.txt")
sbatch --array=0-$((NUM_FILES - 1)) hpc_array_jobs/submit_point_lasso.sh
```

Output: one `.rds` model object per array task, feeding
[`03_combine_results/`](../03_combine_results/).
