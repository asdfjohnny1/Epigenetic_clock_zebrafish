# ==========================================================================
# PIPELINE STAGE: 07 / External validation - step 4: diagnose prediction failure
# ==========================================================================
# The fin clock's predicted-vs-actual correlation on the external dataset
# is close to zero (R^2 = 0.0055) despite good CpG coverage (median 20/22)
# and no chromosome-naming mismatch (both danRer11/UCSC and GRCz11/Ensembl
# are the same assembly and coordinates, just chr-prefixed vs not). This
# re-extracts the external samples' raw methylation values at the clock's
# 22 CpGs and compares their per-CpG mean/SD directly against the training
# data's own per-CpG mean/SD, to check for the classic signature of a
# batch/platform effect: a systematic shift or compression in methylation
# values at the same genomic windows, unrelated to age.
#
# Also tests a second, distinct hypothesis: rather than (or alongside) a
# batch/platform shift, this study's 22 CpGs were selected by elastic net
# because they correlate with age IN THIS STUDY'S OWN COHORT - there is no
# reason to expect the same 22 loci to carry any age signal at all in an
# independent cohort collected for a different clock study (Mayne et al.
# 2020, whose own model uses a different, non-overlapping set of 26-29
# CpGs). This computes, per CpG, its correlation with age within the
# external cohort and compares it against the same CpG's correlation with
# age in the training cohort - directly distinguishing "shifted but still
# age-informative" from "not age-informative in this cohort at all".
#
# Input:  data/Metrics_all_models/best_models/best_model_fin.rds
#         $PROJECT_ROOT/external_validation/cov_files/<sample>.bismark.cov.gz
# Output: results/external_validation/methylation_shift_diagnostic.csv
#         results/external_validation/methylation_shift_diagnostic.png
#         results/external_validation/age_correlation_diagnostic.png
# ==========================================================================

PROJECT_ROOT <- Sys.getenv("EPICLOCK_ROOT", unset = normalizePath("."))

suppressPackageStartupMessages({
  library(here)
  library(methylKit)
  library(GenomicRanges)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(ggrepel)
})

model <- readRDS(here("data/Metrics_all_models/best_models/best_model_fin.rds"))
final_cpg <- model$final_clock$final_cpg

cov_dir <- file.path(PROJECT_ROOT, "external_validation", "cov_files")
meta <- read_csv(here("07_external_validation/external_metadata_template.csv"), show_col_types = FALSE)
output_dir <- here("results/external_validation")

coords <- do.call(rbind, strsplit(gsub("^chr", "", final_cpg), ":"))
window_gr <- GRanges(
  seqnames = coords[, 1],
  ranges   = IRanges(start = as.numeric(coords[, 2]), end = as.numeric(coords[, 3])),
  cpg      = final_cpg
)

extract_methylation_at_windows <- function(cov_file, sample_id, window_gr) {
  meth_raw <- methRead(cov_file, sample.id = sample_id, assembly = "danRer11",
                        treatment = 0, context = "CpG", pipeline = "bismarkCoverage", mincov = 1)
  regional <- regionCounts(meth_raw, window_gr)
  df <- getData(regional)
  full <- tibble(cpg = window_gr$cpg, meth_pct = NA_real_)
  matched <- match(paste(df$chr, df$start, df$end, sep = ":"),
                    paste(as.character(seqnames(window_gr)), start(window_gr), end(window_gr), sep = ":"))
  ok <- !is.na(matched)
  full$meth_pct[matched[ok]] <- ifelse(df$coverage[ok] > 0, 100 * df$numCs[ok] / df$coverage[ok], NA_real_)
  full
}

cov_files <- file.path(cov_dir, paste0(meta$sample_id, ".bismark.cov.gz"))
meta <- meta[file.exists(cov_files), ]
cov_files <- cov_files[file.exists(cov_files)]

message("Re-extracting methylation for ", length(cov_files), " external samples...")
meth_list <- Map(extract_methylation_at_windows, cov_files, meta$sample_id,
                  MoreArgs = list(window_gr = window_gr))
external_mat <- t(sapply(meth_list, function(x) x$meth_pct))
colnames(external_mat) <- final_cpg

train_mat <- model$methylation_data_train.final[, final_cpg, drop = FALSE]

# -----------------------
# Per-CpG mean/SD: training vs. external
# -----------------------
comparison <- tibble(
  cpg           = final_cpg,
  train_mean    = colMeans(train_mat, na.rm = TRUE),
  train_sd      = apply(train_mat, 2, sd, na.rm = TRUE),
  external_mean = colMeans(external_mat, na.rm = TRUE),
  external_sd   = apply(external_mat, 2, sd, na.rm = TRUE),
  external_n    = colSums(!is.na(external_mat))
) %>%
  mutate(mean_shift = external_mean - train_mean)

write_csv(comparison, file.path(output_dir, "methylation_shift_diagnostic.csv"))
message("Wrote: ", file.path(output_dir, "methylation_shift_diagnostic.csv"))
print(comparison %>% select(cpg, train_mean, external_mean, mean_shift) %>% arrange(desc(abs(mean_shift))))

cat("\n=== Overall shift summary ===\n")
cat("Median per-CpG mean methylation, training: ", round(median(comparison$train_mean), 1), "%\n")
cat("Median per-CpG mean methylation, external: ", round(median(comparison$external_mean), 1), "%\n")
cat("Median absolute per-CpG shift:             ", round(median(abs(comparison$mean_shift)), 1), "percentage points\n")
cat("Correlation of per-CpG means (train vs external): ",
    round(cor(comparison$train_mean, comparison$external_mean, use = "complete.obs"), 3),
    "  (near 1 = same relative ranking of CpGs just possibly shifted/compressed;",
    " near 0 = the two datasets don't even agree on which CpGs are hyper/hypomethylated)\n")

# -----------------------
# Plot: per-CpG mean methylation, training vs external, with 1:1 line
# -----------------------
p <- ggplot(comparison, aes(x = train_mean, y = external_mean)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
  geom_point(size = 2, colour = "#C46E7A") +
  geom_text_repel(aes(label = sub("^chr", "", cpg)), size = 2.5, max.overlaps = 30) +
  labs(
    x = "Training-set mean methylation (%)",
    y = "External-dataset mean methylation (%)",
    title = "Per-CpG mean methylation: training vs. external dataset",
    subtitle = "Points on the dashed line = no shift between datasets at that CpG"
  ) +
  coord_equal(xlim = c(0, 100), ylim = c(0, 100)) +
  theme_classic(base_size = 12)

ggsave(file.path(output_dir, "methylation_shift_diagnostic.png"), p, width = 7, height = 7, dpi = 400)
message("Wrote: ", file.path(output_dir, "methylation_shift_diagnostic.png"))

# -----------------------
# Per-CpG correlation with age: training cohort vs. external cohort.
# A CpG can have near-identical mean methylation in both datasets (ruled
# out as a broad batch shift above) while carrying zero age-signal in one
# of them - mean-level agreement and age-informativeness are different
# properties, and only the second one is what the clock actually needs.
# -----------------------
age_col <- grep("^age_(months|weeks|days)$", colnames(meta), value = TRUE)[1]
if (is.na(age_col)) stop("No age_months/age_weeks/age_days column found in metadata.")
age_unit <- sub("^age_", "", age_col)
days_per_month <- 30.44
to_months <- switch(age_unit, months = function(x) x,
                     weeks = function(x) x * 7 / days_per_month,
                     days = function(x) x / days_per_month)
meta <- meta %>% mutate(age_months = to_months(.data[[age_col]]))

train_age <- model$metadata_train.final$Age[match(rownames(train_mat), model$metadata_train.final$Sanger_new_ID)]

age_cor <- tibble(
  cpg = final_cpg,
  train_age_cor    = sapply(seq_along(final_cpg), function(i)
    suppressWarnings(cor(train_mat[, i], train_age, use = "complete.obs"))),
  external_age_cor = sapply(seq_along(final_cpg), function(i)
    suppressWarnings(cor(external_mat[, i], meta$age_months, use = "complete.obs")))
)

comparison <- comparison %>% left_join(age_cor, by = "cpg")
write_csv(comparison, file.path(output_dir, "methylation_shift_diagnostic.csv"))

cat("\n=== Per-CpG correlation with age: training cohort vs. external cohort ===\n")
print(comparison %>% select(cpg, train_age_cor, external_age_cor) %>%
        arrange(desc(abs(train_age_cor))))
cat("\nMedian |correlation with age| in training cohort: ",
    round(median(abs(comparison$train_age_cor), na.rm = TRUE), 2), "\n")
cat("Median |correlation with age| in external cohort: ",
    round(median(abs(comparison$external_age_cor), na.rm = TRUE), 2), "\n")
cat("CpGs with |train_age_cor| > 0.3 that also have |external_age_cor| > 0.3: ",
    sum(abs(comparison$train_age_cor) > 0.3 & abs(comparison$external_age_cor) > 0.3, na.rm = TRUE),
    " of ", sum(abs(comparison$train_age_cor) > 0.3, na.rm = TRUE),
    " CpGs that were actually age-correlated in training\n")

p2 <- comparison %>%
  select(cpg, train_age_cor, external_age_cor) %>%
  pivot_longer(-cpg, names_to = "cohort", values_to = "age_cor") %>%
  mutate(cohort = recode(cohort, train_age_cor = "Training cohort", external_age_cor = "External cohort"),
         cpg = sub("^chr", "", cpg)) %>%
  ggplot(aes(x = reorder(cpg, age_cor), y = age_cor, fill = cohort)) +
  geom_col(position = "dodge") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  coord_flip() +
  labs(x = "", y = "Correlation with chronological age",
       title = "Per-CpG age-correlation: training vs. external cohort", fill = "") +
  theme_classic(base_size = 11) +
  theme(legend.position = "bottom")

ggsave(file.path(output_dir, "age_correlation_diagnostic.png"), p2, width = 8, height = 7, dpi = 400)
message("Wrote: ", file.path(output_dir, "age_correlation_diagnostic.png"))
