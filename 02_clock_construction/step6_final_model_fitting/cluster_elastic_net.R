# ==========================================================================
# PIPELINE STAGE: 02 / Clock construction - step 6 (Cluster, Elastic Net)
# ==========================================================================
# Array-ready: one invocation processes one *_clock_setup.rds scaffold
# (one SLURM array task) produced by the step 1-5 script above.
# Performs the nested 100-outer x 50-inner iterative model fitting
# described in the Methods (candidate CpG subsetting, held-out
# validation, best-outer-iteration selection).
# 
# Submit via the matching script in hpc_array_jobs/.
# 
# Shared verbatim between both aims of the study.
# ==========================================================================

# Set the root directory of this project (edit here, or export EPICLOCK_ROOT
# before launching R / sbatch). All absolute paths below are built from this.
PROJECT_ROOT <- Sys.getenv("EPICLOCK_ROOT", unset = normalizePath("."))

#!/usr/bin/env Rscript
# step4_array_clock.R
# Array-ready Step 6 (Elastic Net): iterative model fitting, testing, and validation.
# Usage:
#   Rscript step4_array_clock.R /path/to/scaffold_file.rds

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(caret)
  library(glmnet)
  library(doParallel)
  library(dplyr)
  library(readxl)
  library(broom)
  library(tibble)     # used minimally; we avoid column_to_rownames to be robust
  library(data.table)
  library(parallel)
})

# --------------------------
# Command line args
# --------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: Rscript step4_array_clock.R <scaffold_file.rds>")

scaffold_file <- args[1]
if (!file.exists(scaffold_file)) stop("Scaffold file does not exist: ", scaffold_file)
message("[", Sys.time(), "] Processing scaffold: ", scaffold_file)

# --------------------------
# Parse path metadata
# --------------------------
# for testing scaffold_file <- file.path(PROJECT_ROOT, "multitissue/making_clocks/R/results/cluster_clock_new/step3.2/Intestines/100/0.2/CClock_methylation_matrix_min20.tsv_step_3_1_step_3_3.rds")


scaffold_dir <- dirname(scaffold_file)
r2_dir_name  <- basename(scaffold_dir)                        
epse_dir    <- dirname(scaffold_dir)                        
epse_name   <- basename(epse_dir )                           
organ_dir     <- basename(dirname(epse_dir))                
organ_name    <- organ_dir
sample_id  <- sub("\\.rds$", "", basename(scaffold_file))

# --------------------------
# Output directories
# --------------------------
output_directory <- file.path(PROJECT_ROOT, "multitissue/making_clocks/R/results/cluster_clock_new")
output_dir4 <- file.path(output_directory, "step4.array.new.new.new/Elastic_net")
metrics_out_dir <- file.path(output_dir4,  organ_name,epse_name, r2_dir_name)
dir.create(metrics_out_dir, recursive = TRUE, showWarnings = FALSE)

# -------------------------------------------------------------------------
# Skip if output for this scaffold already exists
# -------------------------------------------------------------------------
output_file <- file.path(output_dir4, paste0(sample_id, "_clock_modelled.rds"))

if (file.exists(output_file)) {
  message("[", Sys.time(), "] Output already exists: ", output_file)
  message("Skipping scaffold and exiting cleanly.")
  quit(save = "no", status = 0)
}


# --------------------------
# Load metadata
# --------------------------
input_directory_metadata <- file.path(PROJECT_ROOT, "multitissue")
metadata_file <- file.path(input_directory_metadata, "final_metadata_uniqueID_new.xlsx")
if (!file.exists(metadata_file)) {
  warning("Metadata file not found at ", metadata_file, " — continuing but some functionality may fail.")
  sample_data <- data.frame()
} else {
  sample_data <- readxl::read_excel(metadata_file) %>%
    dplyr::select(Clock, sequencing_block, plate_name, Sanger_new_ID, Organ, original_tube_ID, Age, Sex, Treatment,individual_fish_ID) %>%
    dplyr::filter(Clock == "multitissue") %>%
    as.data.frame(stringsAsFactors = FALSE)
}

# --------------------------
# Parallel setup (respect SLURM_CPUS_PER_TASK)
# --------------------------
n_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = parallel::detectCores(logical = TRUE) - 1))
n_cores <- max(1, n_cores)
cl <- makeCluster(n_cores)
doParallel::registerDoParallel(cl)
on.exit({
  try(parallel::stopCluster(cl), silent = TRUE)
}, add = TRUE)
message("Registered parallel backend with ", n_cores, " cores.")

# --------------------------
# Utility: robust predict_age that supports caret glmnet and glmnet/cv.glmnet directly
# model_arg can be caret train object (finalModel) or glmnet / cv.glmnet
# newdata: samples x features (matrix or data.frame)
# --------------------------
predict_age <- function(model_arg, newdata) {
  # If caret::train object (has class 'train'), use predict.train with newdata frame
  if (inherits(model_arg, "train")) {
    preds <- predict(model_arg, newdata = as.data.frame(newdata))
    return(as.numeric(preds))
  }
  # If glmnet or cv.glmnet: predict with newx matrix
  if (inherits(model_arg, "cv.glmnet") || inherits(model_arg, "glmnet")) {
    newx <- as.matrix(newdata)
    # prefer lambda.min if cv.glmnet
    if (inherits(model_arg, "cv.glmnet")) {
      preds <- predict(model_arg, newx = newx, s = "lambda.min")
    } else {
      preds <- predict(model_arg, newx = newx)
    }
    # glmnet predict returns matrix; coerce numeric vector
    preds <- as.numeric(preds)
    return(preds)
  }
  # fallback to generic predict
  preds <- predict(model_arg, newdata = as.data.frame(newdata))
  return(as.numeric(preds))
}

#predict_age <- function(x, y) {
#  predicted_values <- predict(x, newx = y)
# Ensure vector shape (glmnet predict sometimes returns matrix)
# if (is.matrix(predicted_values)) predicted_values <- as.numeric(predicted_values[, 1])
#names(predicted_values) <- NULL
#return(predicted_values)
#}

# Helper: fast column-wise median imputation
median_impute_matrix <- function(mat) {
  mat <- as.matrix(mat)
  na_cols <- colSums(is.na(mat)) > 0
  if (any(na_cols)) {
    medians <- apply(mat, 2, function(x) median(x, na.rm = TRUE))
    for (j in which(na_cols)) {
      mat[is.na(mat[, j]), j] <- medians[j]
    }
  }
  # remove columns entirely NA
  mat <- mat[, colSums(is.na(mat)) < nrow(mat), drop = FALSE]
  return(mat)
}



num_runs_inner <- 100
num_runs_outer <- 50

id_columns <- c("Probe", "chr", "start", "end","strand")


# -----------------------
# Load scaffold and run the same pipeline logic but only for this scaffold file
# -----------------------
process_scaffold <- function(scaffold_file, sample_data, metrics_out_dir, input_directory_base) {
  
  # read scaffold object
  file_entry <- readRDS(scaffold_file)
  
  # Defensive initialisations (in case scaffold missing fields)
  # --- Defensive initialisations
  if (is.null(file_entry$best_metrics)) {
    file_entry$best_metrics <- list(
      best_train_mae = Inf,
      best_test_mae = Inf,
      best_train_r_squared = -Inf,
      best_test_r_squared = -Inf,
      best_val_r_squared = -Inf,
      best_val_mae = Inf,
      best_outter_seed = NA_integer_,
      best_inner_seed = NA_integer_,
      best_nb_cpg = NA_integer_,
      best_cpg_ID = NULL,
      best_probe_list = NULL,
      best_model = NULL
    )
  }
  if (is.null(file_entry$best_models)) file_entry$best_models <- list()
  if (is.null(file_entry$failed_models)) file_entry$failed_models <- list()
  
  # Initialise running mean and SD for online z-scoring ( chosing of best model - aiming for smaller cpg models)
  file_entry$best_metrics$best_score = -Inf 
  r2_mean <- 0
  r2_sd <- 1
  mae_mean <- 0
  mae_sd <- 1
  n <- 0
  cpg_mean <- 0
  cpg_sd <- 1
  beta <- 1 # trade-off parameter for how many cpg are used in mdoel -  a value of 2 tells code clocks with smaller cpg are twice as important as one SD of performance ( R2 and mae) - 1 is equal importance 
  
  
  # If scaffold stores methylation under a different name, try to standardise
  if (is.null(file_entry$methylation_data_clustered) && !is.null(file_entry$methylation_data)) {
    file_entry$methylation_data_clustered <- file_entry$methylation_data
  }
  
  # Populate current_probes
  if (!is.null(file_entry$calculated_metrics) && nrow(file_entry$calculated_metrics) > 0) {
    file_entry$current_probes <- file_entry$calculated_metrics
  } else {
    file_entry$current_probes <- data.frame()
    warning("No calculated_metrics for sample ", sample_id)
  }
  
  # Filter sample metadata to those present in methylation file columns
  sample_cols <- intersect(colnames(file_entry$methylation_data_clustered), sample_data$Sanger_new_ID)
  file_entry$sample_data_filtered <- sample_data %>% dplyr::filter(Sanger_new_ID %in% sample_cols) %>% as.data.frame(stringsAsFactors = FALSE)
  
  
  
  # Main try/catch to protect per-sample processing
  tryCatch({
    
    for (i in 1:num_runs_outer) {
      tryCatch({
        
        message(paste0("Outer loop run ", i, " for R2 threshold ", r2_dir_name))
        set.seed(i)
        outerloop_current_seed <- i
        
        # 6.1: Create training/validation split (80/20)
        if (!is.null(file_entry$sample_data_filtered) &&
            all(c("Age", "Organ") %in% colnames(file_entry$sample_data_filtered))) {
          
          message("Starting STEP 6.1")
          
          meta <- file_entry$sample_data_filtered
          
          # 1. Summarise treatment per individual
          id_summary <- meta |>
            group_by(individual_fish_ID) |>
            summarise(
              primary_treatment = names(sort(table(Treatment), decreasing = TRUE))[1],
              primary_age = names(sort(table(Age), decreasing = TRUE))[1],
              .groups = "drop"
            )
          
          validationIndex <- createDataPartition(
            y = paste(id_summary$primary_treatment, id_summary$primary_age),
            p = 0.8,
            list = FALSE,
            times = 1
          )
          
          train_ids <- id_summary$individual_fish_ID[validationIndex]
          test_ids  <- id_summary$individual_fish_ID[-validationIndex]
          
          # 3. Subset metadata and methylation data
          train_idx <- meta$individual_fish_ID %in% train_ids
          test_idx  <- meta$individual_fish_ID %in% test_ids
          
          
          
          file_entry$validationIndex <- validationIndex
          
          # Final validation metadata (20%)
          file_entry$final_validation_metadata  <- meta[test_idx, , drop = FALSE]
          # Creation clock metadata (80%)
          file_entry$creation_clock_metadata  <- meta[train_idx, , drop = FALSE]
          
          # 4. Sanity check: no leakage and treatment balance
          length(intersect(train_ids, test_ids))  # should be 0
          
          table(file_entry$creation_clock_metadata$Age)
          table(file_entry$final_validation_metadata$Age)
          table(file_entry$creation_clock_metadata$individual_fish_ID)
          table(file_entry$final_validation_metadata$individual_fish_ID)
          
          
          
        } else {
          warning("Missing or incomplete sample_data_filtered for sample ", sample_id)
          next
        }
        
        for (j in 1:num_runs_inner) {
          tryCatch({
            
            message(paste0("Inner loop run ", j, " for R2 threshold ", r2_dir_name))
            set.seed(j)
            innerloop_current_seed <- j
            
            # STEP 6.2 - selecting CPG
            
            
            
            # Randomly select number and subset CpGs
            if (!is.null(file_entry$current_probes) && nrow(file_entry$current_probes) > 1) {
              message("Starting STEP 6.2")
              
              num_candidate_cpg <- sample(1:nrow(file_entry$current_probes), 1)
              candidate_cpg <- file_entry$current_probes[sample(nrow(file_entry$current_probes), num_candidate_cpg), , drop = FALSE]
              
              file_entry$methylation_data_filtered <- file_entry$methylation_data_clustered %>%
                filter(Probe %in% candidate_cpg$Probe)
              
              file_entry$methylation_data_filtered <- column_to_rownames(file_entry$methylation_data_filtered, "Probe")
              file_entry$methylation_data_filtered <- file_entry$methylation_data_filtered %>%
                dplyr::select(-Chromosome, -Start, -End, -Strand, -cluster, -unique_ID)
            } else {
              warning("Missing or incomplete data for = ", organ_name, ",  at r threshold = ", r2_dir_name )
              num_candidate_cpg <- numeric(0)
              candidate_cpg <- numeric(0)
              methylation_data_filtered <- data.frame()
            }
            
            
            
            # 6.3 Clean methylation data
            if (!is.null(file_entry$methylation_data_filtered) &&
                ncol(file_entry$methylation_data_filtered) > 0) {
              message("Starting STEP 6.3")
              
              md_filtered <- file_entry$methylation_data_filtered
              head(rownames(md_filtered))
              md_filtered <- as.matrix(md_filtered)
              # Keep only CpGs on canonical chromosomes (e.g. chr1–chr25)
              canonical_idx <- grepl("^([0-9]+|X|Y|M):", rownames(md_filtered))
              md_filtered <- md_filtered[canonical_idx, , drop = FALSE]
              
              # remove rows whicha re all NA
              md_filtered <- md_filtered[rowSums(is.na(md_filtered)) != ncol(md_filtered), ]
              
              
              # Remove any CpG where more than 20% of the samples (columns) are NA.
              # Compute fraction of NA per CpG 
              #  na_fraction <- rowSums(is.na(md_filtered)) / ncol(md_filtered)
              
              # Keep only rows below threshold
              #  md_filtered <- md_filtered[na_fraction <= 0.1, ]
              
              
              
              # make probe rownames
              #  rownames(md_filtered) <- md_filtered$Probe
              # md_filtered <- md_filtered %>% select(-all_of(id_columns))
              
              #  keep top 10% cpg 
              #  v <- apply(md_filtered, 1, var, na.rm = TRUE)
              # thr <- quantile(v, 0.9)  # keep top 10% most variable CpGs
              # md_filtered <- md_filtered[v > thr, ]
              
              
              file_entry$methylation_data_filtered <- md_filtered
              
            } else {
              warning("Incomplete methylation_data_filtered for sample ", sample_id)
              next
            }
            
            
            
            
            # 6.4 Create final filtered matrix & metadata
            if (!is.null(file_entry$creation_clock_metadata) && nrow(file_entry$creation_clock_metadata) > 0 && nrow(file_entry$methylation_data_filtered) > 0) {
              
              message("Starting STEP 6.4")
              rownames(file_entry$creation_clock_metadata) <- file_entry$creation_clock_metadata$Sanger_new_ID
              
              valid_ids <- intersect(rownames(file_entry$creation_clock_metadata), colnames(file_entry$methylation_data_filtered))
              file_entry$methylation_data_filtered_final <- file_entry$methylation_data_filtered[, valid_ids, drop = FALSE]
              file_entry$methylation_data_filtered_final <- t(file_entry$methylation_data_filtered_final)
              
              file_entry$creation_clock_metadata <- file_entry$creation_clock_metadata[order(rownames(file_entry$creation_clock_metadata)), ]
              file_entry$methylation_data_filtered_final <- file_entry$methylation_data_filtered_final[order(rownames(file_entry$methylation_data_filtered_final)), ]
              
              #
              common_rows <- intersect(rownames(file_entry$creation_clock_metadata), 
                                       rownames(file_entry$methylation_data_filtered_final))
              
              file_entry$creation_clock_metadata <- file_entry$creation_clock_metadata[common_rows, , drop = FALSE]
              file_entry$methylation_data_filtered_final <- file_entry$methylation_data_filtered_final[common_rows, , drop = FALSE]
              
              
              if (!all(rownames(file_entry$creation_clock_metadata) == rownames(file_entry$methylation_data_filtered_final))) {
                stop("Error: Row names are not aligned after sorting!")
              }
            } else {
              warning("Missing or incomplete data for = ", organ_name, ",  at r threshold = ", r2_dir_name )
              file_entry$methylation_data_filtered_final <- data.frame()
            }
            
            
            # 6.5 Split creation_clock_metadata into train/test 75%/25%
            if (!is.null(file_entry$creation_clock_metadata) &&
                ncol(file_entry$methylation_data_filtered_final) > 1 &&
                is.matrix(file_entry$methylation_data_filtered_final)) {
              message("Starting STEP 6.5")
              
              
              meta <- file_entry$creation_clock_metadata
              
              # 1. Summarise treatment per individual
              id_summary <- meta |>
                group_by(individual_fish_ID) |>
                summarise(
                  primary_treatment = names(sort(table(Treatment), decreasing = TRUE))[1],
                  primary_age = names(sort(table(Age), decreasing = TRUE))[1],
                  .groups = "drop"
                )
              
              trainIndex.all <- createDataPartition(
                y = paste(id_summary$primary_treatment, id_summary$primary_age),
                p = 0.75,
                list = FALSE,
                times = 1
              )
              
              train_ids <- id_summary$individual_fish_ID[trainIndex.all]
              test_ids  <- id_summary$individual_fish_ID[-trainIndex.all]
              
              # 3. Subset metadata and methylation data
              train_idx <- meta$individual_fish_ID %in% train_ids
              test_idx  <- meta$individual_fish_ID %in% test_ids
              
              
              # train metadata (75%)
              file_entry$metadata_train.final <- meta[train_idx, , drop = FALSE]
              # test metadata (25%)
              file_entry$metadata_test.final  <- meta[test_idx, , drop = FALSE]  
              
              
              # 4. Sanity check: no leakage and treatment balance
              length(intersect(train_ids, test_ids))  # should be 0
              
              table(file_entry$creation_clock_metadata$Age)
              table(file_entry$final_validation_metadata$Age)
              table(file_entry$creation_clock_metadata$individual_fish_ID)
              table(file_entry$final_validation_metadata$individual_fish_ID)
              
              
              # 5. filter methylation 
              
              
              file_entry$methylation_data_train.final <- file_entry$methylation_data_filtered_final[train_idx, , drop = FALSE]
              file_entry$methylation_data_test.final <- file_entry$methylation_data_filtered_final[test_idx, , drop = FALSE]
              
              if (!all(rownames(file_entry$methylation_data_train.final) == rownames(file_entry$metadata_train.final))) {
                stop("Row names misalignment between train methylation and metadata for sample ", sample_id)
              }
            } else {
              warning("Insufficient data for train/test split for sample ", sample_id)
              next
            }
            
            # 6.6 Remove columns with NA in train and match test columns
            if (nrow(file_entry$methylation_data_train.final) > 0) {
              message("Starting STEP 6.6")
              
              cols_to_keep <- colSums(is.na(file_entry$methylation_data_train.final)) != nrow(file_entry$methylation_data_train.final)
              file_entry$methylation_data_train.final <- file_entry$methylation_data_train.final[, cols_to_keep, drop = FALSE]
              file_entry$methylation_data_test.final <- file_entry$methylation_data_test.final[, colnames(file_entry$methylation_data_test.final) %in% colnames(file_entry$methylation_data_train.final), drop = FALSE]
            } else {
              warning("Empty train data after NA filtering for sample ", sample_id)
              next
            }
            
            # --- STEP 6.7: Elastic Net using caret::train() ---
            if (ncol(file_entry$methylation_data_train.final) > 1) {
              message("STEP 6.7: Elastic Net via caret — sample: ", sample_id)
              
              # Prepare training and test data (ensure same features)
              x_train_full <- median_impute_matrix(file_entry$methylation_data_train.final)
              x_test_full  <- median_impute_matrix(file_entry$methylation_data_test.final)
              common_cols <- intersect(colnames(x_train_full), colnames(x_test_full))
              x_train <- as.matrix(x_train_full[, common_cols, drop = FALSE])
              x_test  <- as.matrix(x_test_full[,  common_cols, drop = FALSE])
              storage.mode(x_train) <- "double"; storage.mode(x_test) <- "double"
              
              # Replace the original
              file_entry$methylation_data_train.final <- x_train
              file_entry$methylation_data_test.final <- x_test
              
              y_train <- as.numeric(file_entry$metadata_train.final$Age)
              y_test  <- as.numeric(file_entry$metadata_test.final$Age)
              
              cv_fit <- cv.glmnet(
                x_train, y_train,
                alpha = 0.5,
                nfolds = 10,
                nlambda = 50,
                type.measure = "mae",
                parallel = TRUE,
                standardize = TRUE
              )
              
              best_alpha <- 0.5
              best_lambda <- cv_fit$lambda.min
              message("Best alpha = ", best_alpha, " | Best lambda = ", signif(best_lambda, 5))
              
              # Predictions
              train_predicted_ages <- predict(cv_fit, newx = x_train, s = best_lambda)
              test_predicted_ages  <- predict(cv_fit, newx = x_test,  s = best_lambda)
              
              train_predicted_ages[train_predicted_ages < 0] <- 0
              test_predicted_ages[test_predicted_ages < 0]   <- 0
              
              # Performance metrics
              train_mae <- mean(abs(y_train - train_predicted_ages))
              test_mae  <- mean(abs(y_test  - test_predicted_ages))
              
              train_r2 <- cor(y_train, train_predicted_ages, use = "complete.obs")^2
              test_r2  <- cor(y_test,  test_predicted_ages,  use = "complete.obs")^2
              
              message("Train MAE = ", round(train_mae,3),
                      " | Test MAE = ", round(test_mae,3),
                      " | Train R² = ", round(train_r2,3),
                      " | Test R² = ", round(test_r2,3))
              
              # Extract coefficients and selected CpGs
              coef_mat <- coef(cv_fit, s = best_lambda)
              coef_vec <- as.numeric(coef_mat)
              coef_names <- rownames(coef_mat)
              selected_idx <- which(coef_vec != 0)
              selected_CpG <- setdiff(coef_names[selected_idx], "(Intercept)")
              n_selected <- length(selected_CpG)
              
              message("Selected CpGs = ", n_selected)
            } else {
              warning("Insufficient CpGs for elastic net model fitting for sample ", sample_id)
            }
            
            
            # 6.8 Keep only best model by z score of both R squared and mae 
            
            if (!is.na(test_r2) &&
                !is.na(test_mae)) {
              message("Starting STEP 6.8")
              
              n <- n+1 
              
              # Update running mean/SD (Welford's algorithm for online variance)
              delta_r2 <- test_r2- r2_mean
              r2_mean <- r2_mean + delta_r2 / n
              delta_mae <- test_mae - mae_mean
              mae_mean <- mae_mean + delta_mae / n
              
              if (n > 1) {
                r2_sd <- sqrt(((n-2)/(n-1)) * r2_sd^2 + delta_r2^2/n)
                mae_sd <- sqrt(((n-2)/(n-1)) * mae_sd^2 + delta_mae^2/n)
              }
              
              epsilon <- 1e-8  # tiny number to avoid division by zero
              
              r2_z  <- (test_r2 - r2_mean) / (r2_sd + epsilon)
              mae_z <- - (test_mae - mae_mean) / (mae_sd + epsilon)
              
              score <- r2_z + mae_z 
              
              if (score >   file_entry$best_metrics$best_score) {
                
                
                file_entry$best_metrics$best_model <- list(
                  score = score,
                  train_mae = train_mae,
                  test_mae = test_mae,
                  train_r2= train_r2,
                  test_r2= test_r2,
                  cv_fit = cv_fit,
                  train_predicted_ages = train_predicted_ages,
                  test_predicted_ages = test_predicted_ages,
                  num_candidate_cpgs =  num_candidate_cpg,
                  candidate_cpg =  candidate_cpg,
                  selected_Cpg =  selected_CpG,
                  best_lambda = best_lambda,
                  best_alpha = best_alpha,
                  metadata_train = file_entry$metadata_train.final,
                  methylation_data_train = file_entry$methylation_data_train.final
                )
                
              } else {
                # Store failed model metrics for tracking
                run_name <- paste0("run_outer_", i, "_inner_", j)
                file_entry$failed_models[[run_name]] <- list(
                  nb_cpg = num_candidate_cpg,
                  cpg_ID = candidate_cpg,
                  selected_Cpg =  selected_Cpg,
                  train_mae = train_mae,
                  test_mae = test_mae,
                  train_r2= train_r2,
                  test_r2= test_r2,
                  probe_ID = r2_dir_name,
                )
              }
            }else {
              warning("R squared & MAE are empty ")
            }
            
          }, error = function(e_inner) {
            message("Inner loop ", j, " failed: ", e_inner$message)
            run_name <- paste0("run_outer_", i, "_inner_", j)
            file_entry$failed_models[[run_name]] <- list(
              error_msg = e_inner$message,
              timestamp = Sys.time()
            )
            # continue to next inner loop iteration
          })
          
        } # end inner loop j
        
        
        # --- STEP 6.9: Refit final clock on full dataset using selected CpGs ---
        
        if (!is.null(file_entry$best_metrics$best_model$cv_fit)) {
          message("Starting STEP 6.9: refitting final elastic net clock on full dataset")
          
          # ---- Extract best lambda and selected CpGs ----
          best_lambda <- file_entry$best_metrics$best_model$best_lambda
          best_alpha  <- file_entry$best_metrics$best_model$best_alpha  # or fixed 0.5
          
          final_cpg <- file_entry$best_metrics$best_model$selected_Cpg
          
          # ---- Combine full dataset ----
          all_data <- rbind(
            file_entry$methylation_data_train.final,
            file_entry$methylation_data_test.final
          )
          all_metadata <- rbind(
            file_entry$metadata_train.final,
            file_entry$metadata_test.final
          )
          
          # ---- Subset and impute ----
          all_data <- all_data[, final_cpg, drop = FALSE]
          na_cols <- colSums(is.na(all_data)) > 0
          if (any(na_cols)) {
            medians <- apply(all_data, 2, median, na.rm = TRUE)
            for (j in which(na_cols)) {
              all_data[is.na(all_data[, j]), j] <- medians[j]
            }
          }
          
          # ---- Fit final condensed clock ----
          final_fit <- glmnet(
            x = as.matrix(all_data),
            y = all_metadata$Age,
            alpha = best_alpha,
            lambda = best_lambda,
            standardize = TRUE
          )
          
          # Predictions
          final_train_predicted_ages <- predict(final_fit, newx =  as.matrix(all_data), s = best_lambda)
          final_train_predicted_ages[final_train_predicted_ages < 0] <- 0
          
          # Performance metrics
          final_train_mae <- mean(abs(all_metadata$Age - final_train_predicted_ages))
          
          final_train_r2 <- cor(all_metadata$Age, final_train_predicted_ages, use = "complete.obs")^2
          
          # ---- Store results ----
          file_entry$final_clock <- list(
            glmnet_fit = final_fit,
            final_cpg = final_cpg,
            n_cpgs = length(final_cpg),
            alpha = best_alpha,
            lambda = best_lambda,
            final_train_predicted_ages = final_train_predicted_ages,
            final_train_mae = final_train_mae,
            final_train_r2 = final_train_r2
          )
          
          message(paste0(
            "Final clock trained with ", length(final_cpg),
            " CpGs (alpha=", round(best_alpha, 2),
            ", lambda=", signif(best_lambda, 3), ")"
          ))
        } else {
          message("No cv.fit mdoel ")
        }
        
        
        
        # STEP 6.10 Validation on the original 20% held out data
        if (!is.null(file_entry$final_clock) && !is.null(file_entry$final_validation_metadata)) {
          message("Starting STEP 6.10")
          
          # Filter validation methylation data
          valid_cols <- intersect(
            colnames(file_entry$methylation_data_clustered),
            file_entry$final_validation_metadata$Sanger_new_ID
          )
          
          
          
          # Subset only valid columns (plus Probe)
          final_validation_methylation_data <- file_entry$methylation_data_clustered %>%
            select(all_of(c("Probe", valid_cols)))
          
          
          
          #  restore probe as row names in a data frame
          final_validation_methylation_data <- as.data.frame(final_validation_methylation_data)
          rownames(final_validation_methylation_data) <- final_validation_methylation_data$Probe
          final_validation_methylation_data$Probe <- NULL
          
          final_validation_methylation_data <- t(final_validation_methylation_data[
            rownames(final_validation_methylation_data) %in% colnames(file_entry$best_metrics$best_model$methylation_data_train),
            , drop = FALSE
          ])
          
          ####### need to make sure validation data is same as final cpg other doesnt work + also need to amke sure columns arent all NA 
          
          # CpGs used in final clock
          final_cpg <- file_entry$best_metrics$best_model$selected_Cpg
          
          # Ensure validation has exactly the same CpGs (in same order)
          valid_mat <- final_validation_methylation_data[, final_cpg, drop = FALSE]
          
          # If some CpGs are missing entirely, add them as NA columns
          missing_cpgs <- setdiff(final_cpg, colnames(valid_mat))
          if (length(missing_cpgs) > 0) {
            for (cpg in missing_cpgs) {
              valid_mat[[cpg]] <- NA
            }
          }
          
          
          
          # Compute medians on training (ignores NA)
          train_meta <- file_entry$best_metrics$best_model$methylation_data_train
          train_meta <- train_meta[, final_cpg, drop = FALSE]
          train_medians <- apply(train_meta, 2, median, na.rm = TRUE)
          
          # Replace NAs in validation with training medians (only if median is not NA)
          for (cpg in names(train_medians)) {
            if (!is.na(train_medians[cpg])) {             # Only use valid medians
              valid_mat[is.na(valid_mat[, cpg]), cpg] <- train_medians[cpg]
            } else {                                      # Optional: fully NA columns get 0
              valid_mat[is.na(valid_mat[, cpg]), cpg] <- 0
            }
          }
          
          final_validation_methylation_data <- valid_mat
          
          ###
          if (ncol(final_validation_methylation_data) > 0) {
            final_val_predicted_ages <- predict( file_entry$final_clock$glmnet_fit, newx = final_validation_methylation_data, s = best_lambda )
            final_validation_mae <- mean(abs(final_val_predicted_ages - file_entry$final_validation_metadata$Age))
            final_validation_r2<- cor(final_val_predicted_ages, file_entry$final_validation_metadata$Age)^2
            
            
          } else {
            warning("No CpGs in validation methylation data for sample ", sample_id)
          }
        }
        
        # 6.11 Keep best validation model
        if (!is.null(file_entry$best_metrics$best_model) && !is.null(final_val_predicted_ages)) {
          
          if (!is.na(final_validation_r2) &&
              final_validation_r2> file_entry$best_metrics$best_val_r_squared) {
            message("Starting STEP 6.11")
            
            file_entry$best_metrics$best_val_r_squared<- final_validation_r2
            file_entry$best_metrics$best_val_mae <- final_validation_mae
            file_entry$best_metrics$final_val_predicted_ages <- final_val_predicted_ages
            file_entry$best_metrics$final_validation_methylation_data <-final_validation_methylation_data
          }
        } else {
          message("Soemthing wrong with best model or validation predictive ages")
        }
        
      }, error = function(e_outer) {
        message("Outer loop ", i, " failed: ", e_outer$message)
        # continue to next outer loop iteration
      })
    } # end outer loop i
    
    # Save updated scaffold with model info
    saveRDS(file_entry, file = file.path(metrics_out_dir, paste0(sample_id, "_clock_modelled.rds")), compress = "xz")
    message("Saved modelled scaffold for sample ", sample_id)
    
  }, error = function(e) {
    message("Error processing sample ", sample_id, ": ", e$message)
    # ---- Save partial scaffold for debugging ----
    partial_file <- sub(
      "step4",
      "step6.array.new.new.new/Elastic_net",
      scaffold_file
    )
    
    # Replace the suffix with "_clock_modelled_error_partial.rds"
    partial_file <- sub(
      "_clock_setup\\.rds$",
      "_clock_modelled_error_partial.rds",
      partial_file
    )
    
    # Ensure the output directory exists
    partial_dir <- dirname(partial_file)
    dir.create(partial_dir, recursive = TRUE, showWarnings = FALSE)
    
    # Try saving the partial object
    try(
      saveRDS(file_entry, file = partial_file, compress = "xz"),
      silent = TRUE
    )
    
    message("Saved partial scaffold for debugging: ", partial_file)
    stop(e)  # rethrow so SLURM sees non-zero exit if desired
  })
}

# Run processing
process_scaffold(scaffold_file, sample_data, metrics_out_dir, input_directory_base)

# Cleanup cluster
parallel::stopCluster(cl)
message("[", Sys.time(), "] Completed: ", sample_id)
