# ==========================================================================
# PIPELINE STAGE: 03 / Combine results (3 of 3)
# ==========================================================================
# Reads the `all_outer_validation_metrics` field written by the step-6
# scripts (02_clock_construction/step6_final_model_fitting/*.R) and
# collates it into one long-format CSV: one row per (model_id, outer
# iteration), instead of the single best-iteration-only row that
# 01_collate_model_metrics.R produces.
#
# This field is only present in *_clock_modelled.rds files produced by the
# current step-6 scripts - re-run step 6 first if your existing .rds
# outputs predate it (check for the field with:
# `!is.null(readRDS(f)$all_outer_validation_metrics)`).
#
# Output consumed by: 05_clock_design_benchmarking/07_check_outer_iteration_
# variability.Rmd, which computes mean +/- 95% CI across the 100 outer
# iterations per model configuration and compares it against the
# single-best-iteration numbers in results_metrics_final.csv.
# ==========================================================================

# Set the root directory of this project (edit here, or export EPICLOCK_ROOT
# before launching R / sbatch). All absolute paths below are built from this.
PROJECT_ROOT <- Sys.getenv("EPICLOCK_ROOT", unset = normalizePath("."))

suppressPackageStartupMessages({
  library(dplyr)
})

input_folders <- list(
  Cluster = file.path(PROJECT_ROOT, "multitissue/making_clocks/R/results/cluster_clock_new/step4.array.new.new.new"),
  Tile    = file.path(PROJECT_ROOT, "multitissue/making_clocks/R/results/tile_clock_final_final/step6.array.new.new"),
  Point   = file.path(PROJECT_ROOT, "multitissue/making_clocks/R/results_split_organ/point_clock/step6.array.new.new/")
)

output_folder <- file.path(PROJECT_ROOT, "multitissue/making_clocks/R/scripts/Graphs/results/final_merged")
dir.create(output_folder, recursive = TRUE, showWarnings = FALSE)
outer_iterations_file <- file.path(output_folder, "results_outer_iterations.csv")

append_or_write <- function(df, filepath) {
  if (!file.exists(filepath)) {
    write.table(df, file = filepath, sep = ",", row.names = FALSE, col.names = TRUE, append = FALSE)
  } else {
    write.table(df, file = filepath, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
  }
}

if (file.exists(outer_iterations_file)) file.remove(outer_iterations_file)

model_counter <- 1
n_missing_field <- 0

for (cpg_grouping in names(input_folders)) {

  input_folder <- input_folders[[cpg_grouping]]

  rds_files <- list.files(
    path = input_folder,
    pattern = "_clock_modelled\\.rds$",
    recursive = TRUE,
    full.names = TRUE
  )

  message("[", cpg_grouping, "] Found ", length(rds_files), " model files")
  pb <- txtProgressBar(min = 0, max = length(rds_files), style = 3)

  for (i in seq_along(rds_files)) {
    file <- rds_files[i]

    tryCatch({
      obj <- readRDS(file)

      if (is.null(obj$all_outer_validation_metrics) || length(obj$all_outer_validation_metrics) == 0) {
        n_missing_field <<- n_missing_field + 1
        setTxtProgressBar(pb, i)
        return(invisible(NULL))
      }

      rel_path <- sub(paste0("^", input_folder, "/"), "", file)
      path_parts <- strsplit(rel_path, "/")[[1]]
      model <- path_parts[1]        # "Lasso" or "Elastic_net"
      organ <- path_parts[2]
      rsq_threshold <- path_parts[4]
      file_stem <- tools::file_path_sans_ext(basename(file))

      filter <- ifelse(
        grepl("min[0-9]+", file_stem),
        sub(".*_(min[0-9]+)(\\.|_).*", "\\1", file_stem),
        "nofilter"
      )

      outer_df <- bind_rows(obj$all_outer_validation_metrics) %>%
        mutate(
          model_id = model_counter,
          cpg_grouping = cpg_grouping,
          model = model,
          organ = organ,
          rsq_threshold = rsq_threshold,
          filter = filter
        ) %>%
        dplyr::select(model_id, cpg_grouping, model, organ, rsq_threshold, filter,
                      outer_iteration, val_r_squared, val_mae, n_cpg)

      append_or_write(outer_df, outer_iterations_file)

      model_counter <- model_counter + 1
      gc()
    }, error = function(e) {
      message("Skipping file: ", file, " - ", e$message)
    })

    setTxtProgressBar(pb, i)
  }
  close(pb)
}

if (n_missing_field > 0) {
  message(n_missing_field, " of the .rds files had no all_outer_validation_metrics field ",
          "(pre-date the step-6 patch) and were skipped.")
}

message("Done. Long-format outer-iteration metrics written to: ", outer_iterations_file)
