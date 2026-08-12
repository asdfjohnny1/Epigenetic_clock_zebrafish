# 05 - Clock design benchmarking

Downstream analysis notebook for Aim 1 of the study: which methodological
choices determine epigenetic clock performance and stability. Corresponds to
the manuscript's Results sections "CpG aggregation strategy...",
"Tissue-specific CpG budgets...", "Feature selection and filtering...", and
"Multi-tissue clocks show decoupled metrics...".

Five numbered parts that must be **run in order within the same R session**
(each part depends on objects created by the previous one - they are not
independently runnable notebooks). A convenience script to knit/run all five
in sequence is left as future work; for now, open them in RStudio in order,
or:

```r
for (f in sort(list.files(pattern = "\\.Rmd$"))) rmarkdown::render(f)
```

| Part | File | Produces |
|---|---|---|
| 1 | `01_setup_load_and_clean.Rmd` | Palettes, helper functions, loads & cleans all 3,036-model result matrices |
| 2 | `02_best_model_ranking_significance.Rmd` | Best model per tissue, significance tests, cross-tissue CpG overlap |
| 3 | `03_determinants_of_performance.Rmd` | LMMs quantifying which design choices drive performance; CpG-count scaling |
| 4 | `04_manuscript_figures_main.Rmd` | **Main figure**: design space & performance, CpG scaling, global design effects (Figure 1 in the manuscript) |
| 5 | `05_manuscript_figures_supplementary.Rmd` | **Supplementary figure**: MAE/R^2 trade-off, preprocessing-threshold robustness, cross-tissue rank stability, tile-window-size analysis (Figure 3 in the manuscript) |

Input: the collated result tables from
[`03_combine_results/`](../03_combine_results/).
