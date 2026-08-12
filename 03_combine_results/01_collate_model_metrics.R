# ==========================================================================
# PIPELINE STAGE: 03 / Combine results (1 of 2)
# ==========================================================================
# Reads every per-model .rds output from step 6 (all tissues x CpG
# groupings x regression families x thresholds) and incrementally writes
# three collated CSVs: metrics (MAE/R^2), predicted ages, and selected
# CpG sets. This is the single results-aggregation step shared by both
# aims of the study - run once, then both aim-specific analyses (05/, 06/)
# read its output.
#
# Output consumed by: 03_combine_results/02_build_summary_tables.R,
# and (indirectly, via that script's output) both aim-specific folders.
# ==========================================================================

# Set the root directory of this project (edit here, or export EPICLOCK_ROOT
# before launching R / sbatch). All absolute paths below are built from this.
PROJECT_ROOT <- Sys.getenv("EPICLOCK_ROOT", unset = normalizePath("."))

# ---------------------------------------------------------------------------
# Script: Plotting MAE and R-squared from epigenetic clocks & collecting model metrics 
## Author: Jean-Charles de Coriolis

# Description:
# This script processes epigenetic clock models saved as .rds files.  
# For each model, it:
#   1. Extracts predicted vs chronological ages for Train, Test, and Validation sets.  
#   2. Collects evaluation metrics (MAE, R²) alongside metadata inferred from the
#      file structure (organ, R² threshold, and CpG filtration level).  
#   3. Generates summary scatterplots of predicted vs chronological age for each set,
#      annotated with performance metrics, and saves them as combined PDFs.  
#   4. Appends the extracted metrics into a single results dataframe, which is exported
#      to CSV for downstream comparative analyses.  
#
# Input:
#   - .rds outputs form STEP 6 of the clock creation code. 
#   -> requires metadata, predicted age, mae & r--squared 
#
# Output:
#   - Combined scatterplots (PDF) for each model.  
#   - A collated CSV file (all_results_metrics.csv) containing MAE and R² values across
#     all organs, thresholds, and filtration settings.  
#
# ---------------------------------------------------------------------------
# ------------------- Setup -------------------
library(dplyr)

input_folder_C  <- file.path(PROJECT_ROOT, "multitissue/making_clocks/R/results/cluster_clock_new/step4.array.new.new.new")
input_folder_T  <- file.path(PROJECT_ROOT, "multitissue/making_clocks/R/results/tile_clock_final_final/step6.array.new.new")
input_folder_P  <- file.path(PROJECT_ROOT, "multitissue/making_clocks/R/results_split_organ/point_clock/step6.array.new.new/")
output_folder <- file.path(PROJECT_ROOT, "multitissue/making_clocks/R/scripts/Graphs/results/final_merged")

okabe_ito <- c(
  "#E69F00", "#56B4E9", "#009E73", "#F0E442",
  "#0072B2", "#D55E00", "#CC79A7", "#999999"
)
okabe_ito_palette <- colorRampPalette(okabe_ito)

# ------------------- Safe write function -------------------
append_or_write <- function(df, filepath) {
  if (!file.exists(filepath)) {
    write.table(df, file = filepath, sep = ",", row.names = FALSE, col.names = TRUE, append = FALSE)
  } else {
    write.table(df, file = filepath, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
  }
}

# ------------------- Collect results -------------------
rds_files <- list.files(
  path = input_folder,
  pattern = "\\.rds$",
  recursive = TRUE,
  full.names = TRUE
)
n_files <- length(rds_files)
safe_value <- function(x) if (is.null(x)) NA else x

metrics_file   <- file.path(output_folder, "cluster_results_metrics.csv")
predicted_file <- file.path(output_folder, "cluster_results_predicted.csv")
cpg_file       <- file.path(output_folder, "cluster_results_cpg.csv")

model_counter <- 1

# --- initialise progress bar ---
pb <- txtProgressBar(min = 0, max = n_files, style = 3)

for (i in seq_along(rds_files)) {
  file <- rds_files[i]
  message("Processing: ", file)
  
  tryCatch({
    obj <- readRDS(file)
    
    rel_path <- sub(paste0("^", input_folder, "/"), "", file)
    path_parts <- strsplit(rel_path, "/")[[1]]
    model <- path_parts[1]
    organ <- path_parts[2]
    epse <- path_parts[3]
    rsq_threshold <- path_parts[4]
    file_stem <- tools::file_path_sans_ext(basename(file))
    tile_methode <- "Cluster"
    
    filter <- ifelse(
      grepl("min[0-9]+", file_stem),
      sub(".*_(min[0-9]+)(\\.|_).*", "\\1", file_stem),
      "nofilter"
    )
    
    # ------------------- Extract lists safely -------------------
    if (model == "Lasso") {
      Train_list <- list(
        data = obj[["best_metrics"]][["best_model"]][["metadata_train"]],
        predicted_ages = obj[["best_metrics"]][["best_model"]][["train_predicted_ages"]],
        mae = obj[["best_metrics"]][["best_model"]][["train_mae"]],
        r_squared = obj[["best_metrics"]][["best_model"]][["train_r_squared"]],
        nb_cpg = obj[["final_clock"]][["n_cpgs"]]
      )
      Validation_list <- list(
        data = obj[["final_validation_metadata"]],
        predicted_ages = obj[["best_metrics"]][["final_val_predicted_ages"]],
        mae = obj[["best_metrics"]][["best_val_mae"]],
        r_squared = obj[["best_metrics"]][["best_val_r_squared"]]
      )
    } else {
      Train_list <- list(
        predicted_ages = obj[["best_metrics"]][["best_model"]][["train_predicted_ages"]],
        mae = obj[["best_metrics"]][["best_model"]][["train_mae"]],
        r_squared = obj[["best_metrics"]][["best_model"]][["train_r2"]],
        nb_cpg = obj[["final_clock"]][["n_cpgs"]]
      )
      Validation_list <- list(
        predicted_ages = obj[["best_metrics"]][["final_val_predicted_ages"]],
        mae = obj[["best_metrics"]][["best_val_mae"]],
        r_squared = obj[["best_metrics"]][["best_val_r_squared"]]
      )
    }
    
    # ------------------- Flatten function -------------------
    flatten_list <- function(lst, dataset_label, model_id) {
      
      metrics_df <- data.frame(
        model_id = model_id,
        dataset = dataset_label,
        model = model,
        organ = organ,
        rsq_threshold = rsq_threshold,
        filter = filter,
        tile_methode = tile_methode,
        epse = epse,
        mae = lst$mae,
        r_squared_value = as.numeric(lst$r_squared),
        stringsAsFactors = FALSE
      )
      
      predicted_df <- cbind(
        data.frame(
          model_id = model_id,
          dataset = dataset_label,
          model = model,
          organ = organ,
          rsq_threshold = rsq_threshold,
          filter = filter,
          tile_methode = tile_methode,
          epse = epse,
          stringsAsFactors = FALSE
        ),
        predicted_age = as.numeric(lst$predicted_ages[, 1]),
        sanger_id = rownames(lst$predicted_ages)
      )
      
      cpg_df <- data.frame(
        model_id = model_id,
        dataset = dataset_label,
        model = model,
        organ = organ,
        rsq_threshold = rsq_threshold,
        filter = filter,
        tile_methode = tile_methode,
        epse = epse,
        num_selected_cpg = Train_list$nb_cpg,
        num_candidate_cpg = obj[["best_metrics"]][["best_model"]][["num_candidate_cpgs"]],
        selected_cpgs = paste(obj[["best_metrics"]][["best_model"]][["selected_Cpg"]], collapse = ";"),
        stringsAsFactors = FALSE
      )
      
      return(list(metrics = metrics_df, predicted = predicted_df, cpg = cpg_df))
    }
    
    # ------------------- Apply to Train & Validation -------------------
    train_flat <- flatten_list(Train_list, "Train", model_counter)
    validation_flat <- flatten_list(Validation_list, "Validation", model_counter)
    
    # ------------------- Stream results to disk -------------------
    append_or_write(train_flat$metrics, metrics_file)
    append_or_write(validation_flat$metrics, metrics_file)
    
    append_or_write(train_flat$predicted, predicted_file)
    append_or_write(validation_flat$predicted, predicted_file)
    
    append_or_write(train_flat$cpg, cpg_file)
    append_or_write(validation_flat$cpg, cpg_file)
    
    model_counter <- model_counter + 1
    
    gc()
  }, error = function(e) {
    message("❌ Skipping file: ", file, " - ", e$message)
  })
  
  # --- update progress bar ---
  setTxtProgressBar(pb, i)
}
close(pb)


message("✅ All results written incrementally to CSV files in: ", output_folder)
