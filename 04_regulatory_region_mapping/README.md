# 04 - Regulatory region mapping

Maps every CpG in the RRBS coverage files to three genomic-context classes
(promoter, gene body, transposable element) for the regulatory-context
analyses used across both aims of the study (methylation density/drift by
genomic context, Polycomb-target enrichment, etc.).

Run in order:

| # | Script | What it does |
|---|---|---|
| 1 | `01_build_regulatory_annotations.sh` | Downloads the Ensembl GRCz11 GTF and derives `promoters_annotated.bed` (TSS +/-1kb) and `gene_bodies_annotated.bed`. |
| 2 | `02_map_cpgs_to_promoters.sh` | SLURM array job: labels each sample's CpGs as promoter-overlapping or not (bedtools). |
| 3 | `03_map_cpgs_to_gene_bodies.sh` | Same, for gene bodies. |
| 4 | `04_map_cpgs_to_transposable_elements.sh` | Builds the shared per-sample file list, then labels CpGs as TE-overlapping or not (requires a RepeatMasker-derived TE BED file - see script header). Run once with no array task set to (re)generate the file list, then submit as an array. |
| 5 | `05_combine_regulatory_region_outputs.sh` | Run once all three array jobs finish. Concatenates every per-sample output into one combined file per region type. Run on an interactive/login node. |
| 6 | `06_regulatory_region_models.Rmd` | Loads the three combined files plus the best-performing clock per tissue (stage 03 output) and produces exploratory regulatory-context figures. |

Scripts 2-5 are standalone, directly runnable SLURM batch/array jobs built
around `bedtools`/`awk`.
