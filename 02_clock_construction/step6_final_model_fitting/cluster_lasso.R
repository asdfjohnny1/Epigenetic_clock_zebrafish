# ==========================================================================
# PIPELINE STAGE: 02 / Clock construction - step 6 (Cluster, Lasso)
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
# Array-ready Step 6 (Lasso): iterative model fitting, testing, and validation.
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
# for testing scaffold_file <- file.path(PROJECT_ROOT, "multitissue/making_clocks/R/results/cluster_clock_new/step3.2/All_organs/1000/0.05/CClock_methylation_matrix_min10.tsv_step_3_1_step_3_3.rds")


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
output_dir4 <- file.path(output_directory, "step4.array.new.new.new/Lasso")
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

predict_age <- function(x, y) {
  predicted_values <- predict(x, newx = y)
  # Ensure vector shape (glmnet predict sometimes returns matrix)
  if (is.matrix(predicted_values)) predicted_values <- as.numeric(predicted_values[, 1])
  names(predicted_values) <- NULL
  return(predicted_values)
}



num_runs_inner <- 100
num_runs_outer <- 50

id_columns <- c("Probe", "chr", "start", "end","strand")


# -----------------------
# Load scaffold and run the same pipeline logic but only for this scaffold file
# -----------------------
process_scaffold <- function(scaffold_file, sample_data, metrics_out_dir) {
  
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
  if (is.null(file_entry$methylation_data_clustered) && !is.null(file_entry$methylation_data_clustered)) {
    file_entry$methylation_data_clustered <- file_entry$methylation_data_clustered
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
              
              
              
              #  keep top 10% cpg 
              #v <- apply(md_filtered, 1, var, na.rm = TRUE)
              # thr <- quantile(v, 0.9, na.rm = TRUE)  # keep top 10% most variable CpGs
              #  md_filtered <- md_filtered[v > thr, ]
              
              
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
            
            # 6.7 Fit lasso model and predict
            if (ncol(file_entry$methylation_data_train.final) > 1) {
              message("Starting STEP 6.7: Fitting lasso model for sample ", sample_id)
              ### we get errors if nas present 
              # so for all_organs nas are replaced with median of feature 
              # -- for train 
              x_train <- file_entry$methylation_data_train.final
              x_train <- as.matrix(x_train)  # ensure it’s a matrix
              
              # Median impute NAs column-wise
              for (i in seq_len(ncol(x_train))) {
                if (anyNA(x_train[, i])) {
                  x_train[, i][is.na(x_train[, i])] <- median(x_train[, i], na.rm = TRUE)
                }
              }
              
              # remove any columsn with just NA 
              x_train <- x_train[, colSums(is.na(x_train)) != nrow(x_train)]
              
              
              
              # -- for test
              x_test <- file_entry$methylation_data_test.final
              x_test <- as.matrix(x_test)  # ensure it’s a matrix
              
              # Median impute NAs column-wise
              for (i in seq_len(ncol(x_test))) {
                if (anyNA(x_test[, i])) {
                  x_test[, i][is.na(x_test[, i])] <- median(x_test[, i], na.rm = TRUE)
                }
              }
              
              # remove any columsn with just NA 
              x_test <- x_test[, colSums(is.na(x_test)) != nrow(x_test)]
              
              # making sure colnames are the same between test and train 
              x_train <- x_train[, colnames(x_train) %in% colnames(x_test)]
              
              
              
              # Replace the original
              file_entry$methylation_data_train.final <- x_train
              file_entry$methylation_data_test.final <- x_test
              
              
              
              # now running model
              cv_fit <- cv.glmnet(
                x = file_entry$methylation_data_train.final,
                y = file_entry$metadata_train.final$Age,
                alpha = 1,
                nfolds = 10,
                parallel = TRUE,
                nlambda = 1000
              )
              
              
              lasso_fit <- glmnet(
                x = file_entry$methylation_data_train.final,
                y = file_entry$metadata_train.final$Age,
                lambda = cv_fit$lambda.min
              )
              
              coefficients <- coef(lasso_fit)
              selected_coeffs <- coefficients[, 1] != 0
              selected_Cpg <-  rownames(coef(lasso_fit))[coef(lasso_fit)[,1] != 0]
              selected_Cpg <- setdiff(selected_Cpg, "(Intercept)")  # remove intercept
              selected_Cpg <- as.character(selected_Cpg)
              
              train_predicted_ages <- predict_age(y = file_entry$methylation_data_train.final, x = lasso_fit)
              test_predicted_ages <- predict_age(y = file_entry$methylation_data_test.final, x = lasso_fit)
              
              train_mae <- mean(abs(train_predicted_ages - file_entry$metadata_train.final$Age))
              test_mae <- mean(abs(test_predicted_ages - file_entry$metadata_test.final$Age))
              
              train_r_squared <- cor(train_predicted_ages, file_entry$metadata_train.final$Age)^2
              test_r_squared <- cor(test_predicted_ages, file_entry$metadata_test.final$Age)^2
              
              
              
              
              
              
            } else {
              warning("Insufficient CpGs for lasso model fitting for sample ", sample_id)
              next
            }
            
            # 6.8 Keep only best model by z score of both R squared and mae 
            
            if (!is.na(test_r_squared) &&
                !is.na(test_mae)) {
              message("Starting STEP 6.8")
              
              n <- n+1 
              
              # Update running mean/SD (Welford's algorithm for online variance)
              delta_r2 <- test_r_squared - r2_mean
              r2_mean <- r2_mean + delta_r2 / n
              delta_mae <- test_mae - mae_mean
              mae_mean <- mae_mean + delta_mae / n
              
              if (n > 1) {
                r2_sd <- sqrt(((n-2)/(n-1)) * r2_sd^2 + delta_r2^2/n)
                mae_sd <- sqrt(((n-2)/(n-1)) * mae_sd^2 + delta_mae^2/n)
              }
              
              # Compute z-scores
              r2_z <- (test_r_squared - r2_mean) / r2_sd
              mae_z <- - (test_mae - mae_mean) / mae_sd  # negative because lower MAE is better
              
              # calculate cpg number score 
              delta_cpg <- num_candidate_cpg - cpg_mean
              cpg_mean <- cpg_mean + delta_cpg / n
              if (n > 1) {
                cpg_sd <- sqrt(((n-2)/(n-1)) * cpg_sd^2 + delta_cpg^2/n)
              }
              
              cpg_z <- - (num_candidate_cpg - cpg_mean) / cpg_sd
              score <- r2_z + mae_z + beta * cpg_z  # equal weighting in z-space while being penalised for smaller cpg number models 
              
              if (score >   file_entry$best_metrics$best_score) {
                
                
                file_entry$best_metrics$best_model <- list(
                  score = score,
                  train_mae = train_mae,
                  test_mae = test_mae,
                  train_r_squared = train_r_squared,
                  test_r_squared = test_r_squared,
                  cv_fit = cv_fit,
                  lasso_fit = lasso_fit,
                  train_predicted_ages = train_predicted_ages,
                  test_predicted_ages = test_predicted_ages,
                  num_candidate_cpgs =  num_candidate_cpg,
                  candidate_cpg =  candidate_cpg,
                  selected_Cpg =  selected_Cpg,
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
                  train_r_squared = train_r_squared,
                  test_r_squared = test_r_squared,
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
        
        if (!is.null(file_entry$best_metrics$best_model)) {
          message("Starting STEP 6.9: refitting final clock on full dataset")
          
          final_cpg <- file_entry$best_metrics$best_model$selected_Cpg
          all_data <- rbind(
            file_entry$methylation_data_train.final,
            file_entry$methylation_data_test.final
          )
          all_metadata <- rbind(
            file_entry$metadata_train.final,
            file_entry$metadata_test.final
          )
          
          # Subset to selected CpGs
          all_data <- all_data[, final_cpg, drop = FALSE]
          
          # Median imputation again if needed
          for (i in seq_len(ncol(all_data))) {
            if (anyNA(all_data[, i])) {
              all_data[, i][is.na(all_data[, i])] <- median(all_data[, i], na.rm = TRUE)
            }
          }
          
          # Fit final model
          final_fit <- glmnet(
            x = as.matrix(all_data),
            y = all_metadata$Age,
            alpha = 1,
            lambda = file_entry$best_metrics$best_model$cv_fit$lambda.min
          )
          
          file_entry$final_clock <- list(
            glmnet_fit = final_fit,
            final_cpg = final_cpg,
            n_cpgs = length(final_cpg),
            lambda = file_entry$best_metrics$best_model$cv_fit$lambda.min
          )
          
          
        }
        
        # STEP 6.10 Validation on the original 20% held out data
        if (!is.null(file_entry$final_clock) && !is.null(file_entry$final_validation_metadata)) {
          message("Starting STEP 6.10")
          
          
          # Identify the valid columns
          valid_cols <- intersect(
            colnames(file_entry$methylation_data_clustered),
            file_entry$final_validation_metadata$Sanger_new_ID
          )
          
          # Subset only 'Probe' and the valid columns
          final_validation_methylation_data <- file_entry$methylation_data_clustered[, c("Probe", valid_cols), with = FALSE]
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
            final_val_predicted_ages <- predict_age(y = final_validation_methylation_data, x =  file_entry$final_clock$glmnet_fit)
            
            final_validation_mae <- mean(abs(final_val_predicted_ages - file_entry$final_validation_metadata$Age))
            final_validation_r_squared <- cor(final_val_predicted_ages, file_entry$final_validation_metadata$Age)^2
            
            
          } else {
            warning("No CpGs in validation methylation data for sample ", sample_id)
          }
        }
        
        # 6.11 Keep best validation model
        if (!is.null(file_entry$best_metrics$best_model) && !is.null(final_val_predicted_ages)) {
          
          if (!is.na(final_validation_r_squared) &&
              final_validation_r_squared > file_entry$best_metrics$best_val_r_squared) {
            message("Starting STEP 6.11")
            
            file_entry$best_metrics$best_val_r_squared <- final_validation_r_squared
            file_entry$best_metrics$best_val_mae <- final_validation_mae
            file_entry$best_metrics$final_val_predicted_ages <- final_val_predicted_ages
            file_entry$best_metrics$final_validation_methylation_data <-final_validation_methylation_data
          }
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
      "step3.2",
      "step4.array.new.new.new/Lasso",
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
process_scaffold(scaffold_file, sample_data, metrics_out_dir)

# Cleanup cluster
parallel::stopCluster(cl)
message("[", Sys.time(), "] Completed: ", sample_id)
