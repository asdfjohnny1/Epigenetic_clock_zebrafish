# ==========================================================================
# PIPELINE STAGE: 07 / External validation - step 3: apply trained clocks
# ==========================================================================
# Applies each tissue's trained final clock (data/Metrics_all_models/
# best_models/best_model_*.rds) to a new, independent dataset - extracting
# methylation at exactly the CpGs/windows each clock was trained on from
# the external dataset's Bismark .cov.gz coverage files (produced by
# 02_bismark_align_extract.sh), predicting age, and, if chronological ages
# are available for the external samples, reporting validation MAE/R^2
# alongside this study's own validation performance for direct comparison.
#
# Input:  data/Metrics_all_models/best_models/best_model_<tissue>.rds
#         $PROJECT_ROOT/external_validation/cov_files/<sample>.bismark.cov.gz
#         07_external_validation/external_metadata_template.csv (filled in
#           with real sample IDs; the age column is optional - leave blank
#           for samples with unknown chronological age)
# Output: results/external_validation/predicted_ages.csv
#         results/external_validation/validation_summary.csv (only for
#           tissues where at least one sample has a known age)
#         results/external_validation/predicted_vs_actual.png
#
# UNITS: this study's own Age metadata (and therefore every trained clock's
# predictions) is in MONTHS - confirmed against data/final_metadata_
# uniqueID_new.xlsx, where the "multitissue" cohort's Age column ranges
# 3-22. glmnet predicts on whatever scale its training y was in, so
# predicted_age below comes out in months regardless of the external
# dataset's own units. external_metadata_template.csv has an age_unit
# column for exactly this reason (Mayne et al.'s public dataset reports
# age in weeks, not months) - ages are converted to months below before
# any comparison is made, so MAE/R^2 and the predicted-vs-actual plot are
# both on a consistent scale. Getting this wrong doesn't silently break
# R^2 (Pearson correlation is scale-invariant) but does inflate/distort
# MAE and makes the plot's 1:1 reference line meaningless.
# ==========================================================================

# Set the root directory of this project (edit here, or export EPICLOCK_ROOT
# before running).
PROJECT_ROOT <- Sys.getenv("EPICLOCK_ROOT", unset = normalizePath("."))

suppressPackageStartupMessages({
  library(here)
  library(methylKit)
  library(GenomicRanges)
  library(glmnet)
  library(dplyr)
  library(readr)
  library(ggplot2)
})

# -----------------------
# Configuration
# -----------------------
# Which tissue(s) to validate against - match to whichever tissue(s) the
# external dataset actually samples. best_model_all_organs.rds is excluded
# from this repository (see data/README.md) - regenerate it via
# 03_combine_results/ if you need the Multi_organ clock, or drop it below.
best_model_files <- list(
  Brain       = here("data/Metrics_all_models/best_models/best_model_brain.rds"),
  Caudal_Fin  = here("data/Metrics_all_models/best_models/best_model_fin.rds"),
  Intestines  = here("data/Metrics_all_models/best_models/best_model_intestines.rds"),
  Multi_organ = here("data/Metrics_all_models/best_models/best_model_all_organs.rds")
)
best_model_files <- best_model_files[file.exists(unlist(best_model_files))]
stopifnot(length(best_model_files) > 0)

cov_dir <- file.path(PROJECT_ROOT, "external_validation", "cov_files")
metadata_file <- here("07_external_validation/external_metadata_template.csv")

output_dir <- here("results/external_validation")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------
# Load external sample metadata
# -----------------------
# Expected columns: sample_id, tissue, plus exactly one age column named
# age_months / age_weeks / age_days (unit is read from the column name
# itself, so it can't be forgotten or set inconsistently) - blank/NA age
# cells still get a predicted age, just no validation metric. Any other
# columns (e.g. sex) are ignored here but pass through harmlessly.
meta <- read_csv(metadata_file, show_col_types = FALSE)
stopifnot(all(c("sample_id", "tissue") %in% colnames(meta)))

age_col <- grep("^age_(months|weeks|days)$", colnames(meta), value = TRUE)
if (length(age_col) != 1) {
  stop("Expected exactly one column named age_months / age_weeks / age_days ",
       "in ", metadata_file, ", found: ",
       if (length(age_col) == 0) "none" else paste(age_col, collapse = ", "))
}
age_unit <- sub("^age_", "", age_col)
message("Reading chronological age from column '", age_col, "' (unit: ", age_unit, ")")

days_per_month <- 30.44   # 365.25 / 12
to_months <- switch(age_unit,
  months = function(x) x,
  weeks  = function(x) x * 7 / days_per_month,
  days   = function(x) x / days_per_month
)
meta <- meta %>% mutate(age_months = to_months(.data[[age_col]]))

# -----------------------
# Helper: read one sample's coverage file restricted to a set of CpG/tile
# windows, returning percent methylation per window (NA if uncovered)
# -----------------------
extract_methylation_at_windows <- function(cov_file, sample_id, window_gr) {
  meth_raw <- methRead(
    location  = cov_file,
    sample.id = sample_id,
    assembly  = "danRer11",
    treatment = 0,
    context   = "CpG",
    pipeline  = "bismarkCoverage",
    mincov    = 1
  )

  regional <- regionCounts(meth_raw, window_gr)
  df <- getData(regional)

  # regionCounts drops windows with zero coverage in this sample - reindex
  # onto the full window set so every clock feature gets a column, NA where
  # uncovered (matches window_gr's own "cpg" mcol id, e.g. "chr3:100:200")
  full <- tibble(cpg = window_gr$cpg, meth_pct = NA_real_)
  matched <- match(paste(df$chr, df$start, df$end, sep = ":"),
                    paste(as.character(seqnames(window_gr)), start(window_gr), end(window_gr), sep = ":"))
  ok <- !is.na(matched)
  full$meth_pct[matched[ok]] <- ifelse(df$coverage[ok] > 0,
                                        100 * df$numCs[ok] / df$coverage[ok], NA_real_)
  full
}

# -----------------------
# Apply one tissue's clock to all available external samples
# -----------------------
apply_one_clock <- function(tissue, model_file) {

  message("=== ", tissue, " ===")
  model <- readRDS(model_file)
  glmnet_fit <- model$final_clock$glmnet_fit
  lambda     <- model$final_clock$lambda
  final_cpg  <- model$final_clock$final_cpg   # e.g. "chr3:5080001:5082000"

  # Coordinates -> GRanges. final_cpg carries a "chr" prefix; Bismark
  # coverage files from a fresh alignment may or may not, depending on the
  # genome FASTA's header naming - normalise both sides to no-prefix here
  # and re-add "chr" only for the mcol id used to match back against them.
  coords <- do.call(rbind, strsplit(gsub("^chr", "", final_cpg), ":"))
  window_gr <- GRanges(
    seqnames = coords[, 1],
    ranges   = IRanges(start = as.numeric(coords[, 2]), end = as.numeric(coords[, 3])),
    cpg      = final_cpg
  )

  samples_this_tissue <- meta %>% filter(tissue == !!tissue)
  if (nrow(samples_this_tissue) == 0) {
    message("No external samples listed for ", tissue, " - skipping.")
    return(NULL)
  }

  cov_files <- file.path(cov_dir, paste0(samples_this_tissue$sample_id, ".bismark.cov.gz"))
  missing <- !file.exists(cov_files)
  if (any(missing)) {
    warning(sum(missing), " coverage file(s) not found for ", tissue,
            " - check 02_bismark_align_extract.sh ran for: ",
            paste(samples_this_tissue$sample_id[missing], collapse = ", "))
  }
  samples_this_tissue <- samples_this_tissue[!missing, ]
  cov_files <- cov_files[!missing]
  if (nrow(samples_this_tissue) == 0) return(NULL)

  meth_list <- Map(extract_methylation_at_windows,
                    cov_files, samples_this_tissue$sample_id,
                    MoreArgs = list(window_gr = window_gr))

  meth_mat <- sapply(meth_list, function(x) x$meth_pct)
  rownames(meth_mat) <- final_cpg
  meth_mat <- t(meth_mat)  # samples x CpGs, matching training orientation

  # Median-impute missing values using this study's own training-set
  # medians (mirrors STEP 6.10 in 02_clock_construction/step6_final_model_
  # fitting/*.R, so uncovered external CpGs are handled the same way
  # held-out validation samples are handled internally)
  train_meth <- model$methylation_data_train.final[, final_cpg, drop = FALSE]
  train_medians <- apply(train_meth, 2, median, na.rm = TRUE)
  for (cpg in colnames(meth_mat)) {
    na_idx <- is.na(meth_mat[, cpg])
    if (any(na_idx)) {
      meth_mat[na_idx, cpg] <- if (!is.na(train_medians[cpg])) train_medians[cpg] else 0
    }
  }

  predicted_age_months <- as.numeric(predict(glmnet_fit, newx = meth_mat, s = lambda))
  predicted_age_months[predicted_age_months < 0] <- 0

  tibble(
    sample_id            = samples_this_tissue$sample_id,
    tissue               = tissue,
    !!age_col            := samples_this_tissue[[age_col]],
    actual_age_months    = samples_this_tissue$age_months,
    predicted_age_months = predicted_age_months,
    n_cpg_used           = model$final_clock$n_cpgs,
    n_cpg_covered        = rowSums(!is.na(t(sapply(meth_list, function(x) x$meth_pct))))
  )
}

results_all <- bind_rows(lapply(names(best_model_files), function(t) {
  apply_one_clock(t, best_model_files[[t]])
}))

write_csv(results_all, file.path(output_dir, "predicted_ages.csv"))
message("Wrote: ", file.path(output_dir, "predicted_ages.csv"))

# -----------------------
# Validation metrics, where actual age is known - both ages compared in
# months (this study's native unit; see the UNITS note at the top of this
# script) against this study's own reported validation performance
# (best_metrics$best_val_mae / best_val_r_squared, already in each
# best_model_*.rds)
# -----------------------
validated <- results_all %>% filter(!is.na(actual_age_months))

if (nrow(validated) > 0) {

  external_summary <- validated %>%
    group_by(tissue) %>%
    summarise(
      n_samples          = n(),
      external_mae       = mean(abs(actual_age_months - predicted_age_months)),
      external_r_squared = suppressWarnings(cor(actual_age_months, predicted_age_months, use = "complete.obs")^2),
      .groups = "drop"
    )

  own_val_metrics <- bind_rows(lapply(names(best_model_files), function(t) {
    m <- readRDS(best_model_files[[t]])
    tibble(tissue = t,
           own_val_mae       = m$best_metrics$best_val_mae,
           own_val_r_squared = m$best_metrics$best_val_r_squared)
  }))

  comparison <- external_summary %>% left_join(own_val_metrics, by = "tissue")
  write_csv(comparison, file.path(output_dir, "validation_summary.csv"))
  message("Wrote: ", file.path(output_dir, "validation_summary.csv"),
          " (external_mae/external_r_squared are both computed in months - ",
          "see the UNITS note at the top of this script)")
  print(comparison)

  p <- ggplot(validated, aes(x = actual_age_months, y = predicted_age_months, colour = tissue)) +
    geom_point(alpha = 0.6, size = 2) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    geom_smooth(method = "lm", se = FALSE) +
    facet_wrap(~tissue) +
    labs(x = "Chronological age (months)", y = "Predicted age (months)",
         title = "External validation: predicted vs. actual age") +
    theme_classic(base_size = 12) +
    theme(legend.position = "none")

  ggsave(file.path(output_dir, "predicted_vs_actual.png"), p, width = 10, height = 5, dpi = 400)
  message("Wrote: ", file.path(output_dir, "predicted_vs_actual.png"))

} else {
  message("No external samples have a known chronological age (empty 'age' column in ",
          "external_metadata_template.csv) - predicted ages were still written to ",
          "predicted_ages.csv, but no validation MAE/R^2 could be computed.")
}
