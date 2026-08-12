# Epigenetic clocks and tissue-specific ageing in zebrafish

Code accompanying a study of DNA-methylation-based epigenetic clocks in
zebrafish (*Danio rerio*), addressing two complementary aims:

1. **Clock design benchmarking** - which methodological choices (CpG
   aggregation strategy, feature filtering, regression family, tissue scope)
   determine epigenetic clock performance and stability.
2. **Cross-tissue ageing biology** - what the best-performing clocks reveal
   about tissue-specific vs. shared ("mosaic") epigenetic ageing across
   brain, caudal fin and intestine.

Both aims draw on the same underlying dataset: 182 RRBS samples (62 brain,
61 caudal fin, 59 intestine) and the same 3,036 trained epigenetic clock
models. The pipeline below reflects that: stages 01-04 are shared
data-processing and clock-construction steps common to both aims, and
stages 05-06 are their respective downstream analyses.

> Manuscript status: in preparation for submission to Nature Aging. See
> `CITATION.cff` for the current working title; update the citation and
> repository URL once submitted/published.

## Pipeline overview

| Stage | Folder | What it does |
|---|---|---|
| 1 | [`01_annotation/`](01_annotation/) | Build per-CpG methylation matrices from Bismark coverage files |
| 2 | [`02_clock_construction/`](02_clock_construction/) | Train all 4,320 candidate clock configurations (steps 1-6) |
| 3 | [`03_combine_results/`](03_combine_results/) | Collate per-model outputs into master result tables (3,036 models pass QC) |
| 4 | [`04_regulatory_region_mapping/`](04_regulatory_region_mapping/) | Map CpGs to promoters / gene bodies / transposable elements |
| 5 | [`05_clock_design_benchmarking/`](05_clock_design_benchmarking/) | Aim 1 downstream analysis & figures |
| 6 | [`06_cross_tissue_ageing_biology/`](06_cross_tissue_ageing_biology/) | Aim 2 downstream analysis & figures |

Each folder has its own README with a script-by-script breakdown. Run stages
in numeric order; within stages 05/06, run the numbered `.Rmd` files in order
within one R session (they build on each other's environment rather than
each reloading everything from scratch - see those folders' READMEs).

Stages 01-04 are shared between both aims.

## Getting started

```bash
git clone <this-repo>
cd epigenetic-clock-zebrafish

# 1. Install dependencies - see environment/R_dependencies.txt
Rscript -e 'source("environment/R_dependencies.txt"); install.packages(cran_packages); BiocManager::install(bioc_packages)'

# 2. Point everything at your own working directory
export EPICLOCK_ROOT=/path/to/your/project/root

# 3. Run stages in order (adjust #SBATCH directives for your own cluster)
```

Most small/derived datasets needed for stages 03 onward (metadata, collated
model metrics, best-model objects, etc.) are bundled in `data/`. Raw
sequencing data and a handful of large intermediate files are not - see
[`data/README.md`](data/README.md) for the full inventory, the ENA
accession, and exactly what won't run without them.

## Environment

R 4.5.x; full CRAN/Bioconductor package list in
[`environment/R_dependencies.txt`](environment/R_dependencies.txt). Shell
pipeline stages additionally require `bedtools`, `gawk`, Trim Galore and
Bismark (see that file for details). All `*.sh` scripts are written as SLURM
batch/array jobs - adapt the `#SBATCH` directives (or the submission method
entirely) for your own scheduler.

## Repository conventions

- Every script/notebook opens with a `PIPELINE STAGE` banner comment stating
  what it does and what it depends on having already been run.
- All absolute paths are built from a `PROJECT_ROOT`/`EPICLOCK_ROOT`
  variable - set `EPICLOCK_ROOT` (or edit the default at the top of each
  script) to your own path.
- The aim-specific analysis notebooks in stages 05/06 are split into
  numbered, topic-scoped `.Rmd` files, meant to be run in order within one R
  session (see each folder's README for the run order).

## Data availability

Raw RRBS sequencing data: ENA accession **PRJEB107117**
(https://www.ebi.ac.uk/ena/browser/view/PRJEB107117), under embargo until
acceptance for publication. Derived data needed to reproduce the figures
without re-running the full pipeline is bundled in [`data/`](data/README.md).

## Citation

See [`CITATION.cff`](CITATION.cff).

## License

MIT - see [`LICENSE`](LICENSE).

---

