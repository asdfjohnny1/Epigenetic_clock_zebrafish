# 06 - Cross-tissue ageing biology

Downstream analysis notebook for Aim 2 of the study: what the
best-performing clocks reveal about tissue-specific vs. shared epigenetic
ageing. Corresponds to the manuscript's Results sections "Tissue-specific
epigenetic age acceleration...", "Tissue-specific and shared CpG
architectures...", and "Methylation density, stochastic drift and entropy...".

Ten numbered parts that must be **run in order within the same R session**
(each part depends on objects created by the previous one). See
[`05_clock_design_benchmarking/README.md`](../05_clock_design_benchmarking/README.md)
for the render-in-sequence snippet.

| Part | File | Produces |
|---|---|---|
| 1 | `01_setup_palettes_functions_data_loading.Rmd` | Palettes, plot theme, best-model and methylation data loading, helper functions |
| 2 | `02_quality_control_and_model_matrices.Rmd` | Batch-effect PCA/QC; loads & cleans all 3,036-model result matrices |
| 3 | `03_best_model_ranking_significance.Rmd` | Best model per tissue, significance tests, cross-tissue CpG overlap |
| 4 | `04_delta_age_trajectories.Rmd` | Epigenetic age acceleration (Delta-Age) by tissue and across chronological age |
| 5 | `05_cpg_methylation_dynamics.Rmd` | CpG methylation trajectories (age/sex) for tissue-specific and overlapping CpG sets; heatmap statistics |
| 6 | `06_cpg_functional_annotation_enrichment.Rmd` | GO/KEGG enrichment and genomic context of tissue-specific and universal overlapping CpGs |
| 7 | `07_methylation_drift_and_pcgt_enrichment.Rmd` | Directionality of age-related methylation drift; Polycomb Group Target (PCGT) gene enrichment |
| 8 | `08_dynamic_network_pca_entropy.Rmd` | Tile-/tissue-level PCA; methylation variance and entropy per tissue |
| 9 | `09_regulatory_regions_analysis.Rmd` | Regulatory-context (promoter/gene-body/TE) methylation summary; uses output of [`04_regulatory_region_mapping/`](../04_regulatory_region_mapping/) |
| 10 | `10_manuscript_figures.Rmd` | **Main figures 2, 4, 5, 6, 7** in the manuscript (clock performance by tissue/CpG grouping; Delta-Age; tissue-specific CpG methylation; overlapping-CpG methylation; system-level divergence) plus supplementary figures 1-3 |

Input: the collated result tables from
[`03_combine_results/`](../03_combine_results/), and the mapped regulatory
regions from [`04_regulatory_region_mapping/`](../04_regulatory_region_mapping/)
(needed for part 9 only).

## Supplementary checks

Three standalone, self-contained scripts (not part of the numbered 1-10
sequence above, and not dependent on it) that re-check specific figure/text
claims. All three are fully runnable from data already tracked in this repo
(`results_cpg_final.csv`, `methylation_cpgs.csv`) - none needs the excluded
large per-organ methylation matrices.

| File | Checks |
|---|---|
| `11_check_drift_duplication_and_direction.Rmd` | Whether Figures 4D and 5D are drawn from the same underlying CpG set, and re-derives the per-organ drift direction for all three CpG-set definitions used across the paper. |
| `12_check_entropy_direction.Rmd` | Re-derives the CpG-wise and sample-wise entropy ranking per organ behind Figure 7C. |
| `13_check_exact_panel_ageing_rate.Rmd` | Re-runs the clock-free methylation-slope ageing-rate check (Kruskal-Wallis + pairwise Wilcoxon) on the exact 30-CpG overlapping panel used elsewhere in the paper, instead of the broader 170-CpG "robust CpG" panel. |
