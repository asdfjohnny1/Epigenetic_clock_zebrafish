# 01 - Annotation

Builds the per-CpG methylation matrices used by every downstream clock
(Point + Tile CpG-island variants) from Bismark `.cov.gz` coverage files.

| Script | What it does |
|---|---|
| `making_annotation.R` | Reads Bismark coverage files, applies coverage filtering, unites samples (CpG present in >=90% of samples), and writes methylation-percentage matrices at multiple filtering stringencies and Tile window sizes. |
| `submit_making_annotation.sh` | SLURM submission wrapper for the above (single large-memory job, not an array). |

**Run first.** Output feeds [`02_clock_construction/`](../02_clock_construction/).

Shared identically between both aims of the study (clock design benchmarking
and cross-tissue ageing biology) - both use the exact same underlying
methylation matrices.

```bash
export EPICLOCK_ROOT=/path/to/your/project/root
sbatch submit_making_annotation.sh
```
