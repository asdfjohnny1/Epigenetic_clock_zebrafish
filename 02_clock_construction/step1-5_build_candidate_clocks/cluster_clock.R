# ==========================================================================
# PIPELINE STAGE: 02 / Clock construction - steps 1-5 (Cluster)
# ==========================================================================
# Builds candidate Cluster epigenetic clocks: data prep, probe filtering,
# per-probe R^2 computation, and modelling-scaffold construction per
# tissue x R^2-threshold combination. Reads annotation output from stage 01.
# 
# Step 6 (final iterative model fitting) is a separate array job -
# see 02_clock_construction/step6_final_model_fitting/.
# 
# Shared verbatim between both aims of the study.
# ==========================================================================

# Set the root directory of this project (edit here, or export EPICLOCK_ROOT
# before launching R / sbatch). All absolute paths below are built from this.
PROJECT_ROOT <- Sys.getenv("EPICLOCK_ROOT", unset = normalizePath("."))

# ========================================================================
# Title: Elastic Net Epigenetic Clock Construction (Organ-wise, Parallelised)
# Author: Jean-Charles de Coriolis

# Description:
# This script performs elastic net regression-based epigenetic clock training 
# for multiple organs using methylation annotation matrices. It is optimised 
# for high-performance computing (HPC) environments and supports parallelised 
# model training and memory-efficient output handling.
#
# Key Features:
# - Loads all required packages with checks
# - Sets up parallel execution using available cores
# - Iterates through input directories corresponding to different organs
# - Trains `num_runs` elastic net models per organ using `caret::train`
# - Saves each trained model to disk immediately to avoid memory overload
# - Designed for robust large-scale epigenetic clock construction
#
# Expected Input:
# - For each organ: an RDS file containing a methylation matrix with sample 
#   ages and CpG/tile features
#
# Output:
# - Per-organ directories containing RDS files of trained models (one per run)
#
# Dependencies:
# - caret, glmnet, doParallel, glmnetUtils, tidyverse, and others
#
# Usage:
# - Configure `input_directory_base`, `output_directory`, and `num_runs`
# - Adapt model formula and data loading path (`your_input_file.rds`) as needed
#
# ========================================================================


# -----------------------------------------------------
# Load Libraries
# -----------------------------------------------------
required_packages <- c(
  "caret", "glmnet", "doMC", "glmnetUtils", "tidyverse", "readxl",
  "GenomicRanges", "data.table", "dplyr", "tidyr", "stringr", "dbscan",
  "progress", "doParallel", "shades", "RColorBrewer", "gridExtra",
  "patchwork", "progressr", "purrr"
)
install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cran.rstudio.com/")
  }
}
invisible(lapply(required_packages, install_if_missing))

invisible(lapply(required_packages, function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf("Package '%s' is not installed. Please install it before running.", pkg))
  }
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}))

# -----------------------------------------------------
# Parallel Setup
# -----------------------------------------------------
n_cores <- parallel::detectCores() - 1
cl <- makeCluster(n_cores)
registerDoParallel(cl)

# Define variables (user must set these accordingly)

# List organ folders


output_directory <- file.path(PROJECT_ROOT, "multitissue/making_clocks/R/results/cluster_clock_new")
input_directory_base <- file.path(PROJECT_ROOT, "multitissue/making_annotation_files/R/Results/cluster_clock")
eps_cluster_dir <- file.path(output_directory, "EPS_data")
input_directory_metadata <-  file.path(PROJECT_ROOT, "multitissue")

# List organ folders
organ_dirs <- list.dirs(path = input_directory_base, recursive = FALSE, full.names = TRUE)
organ_names <- basename(organ_dirs)

# laod in metadata 
sample_data <- read_excel(file.path(input_directory_metadata, "final_metadata_uniqueID_new.xlsx")) %>%
  select(Clock,sequencing_block, plate_name, Sanger_new_ID, Organ, original_tube_ID, Age, Sex, Treatment) %>%
  filter(Clock == "multitissue")


# Define ID columns
id_columns <- c("Probe", "Chromosome", "Start", "End", "Strand")


num_runs <- 100

# Set input/output directories




# STEP 2: Create file list for each organ
create_file_list <- function(organ_folder_path, organ_name) {
  message("Processing organ: ", organ_name)
  
  file.list <- list.files(
    path = organ_folder_path,
    pattern = "*.tsv",
    full.names = TRUE
  )
  
  if (length(file.list) == 0) {
    warning("No matching files found in ", organ_folder_path)
    return(NULL)
  }
  
  return(file.list)
}

file_list <- setNames(
  lapply(seq_along(organ_dirs), function(i) {
    dir_path <- organ_dirs[i]
    org_name <- organ_names[i]
    
    if (!dir.exists(dir_path)) {
      warning("Skipping missing directory: ", dir_path)
      return(NULL)
    }
    
    create_file_list(dir_path, org_name)
  }),
  organ_names
)

# Step 3: Define parameters for clustering and filtering
epses <- c("22500", "3000", "2000", "1500", "1000", "750", "500", "100", "50")
thresholds <- c("nofilter", "0.05", "0.1", "0.2", "0.3")
organ_names <- names(file_list)

# Predict age function
predict_age <- function(x, y) {
  predicted_values <- predict(x, newx = y)
  colnames(predicted_values) <- "Predicted_Age"
  return(predicted_values)
}


# for testing 
#organ <- organ_names[1]
#organ_specific_file_list <- file_list[[organ]]  
#file_path <- organ_specific_file_list[1]
#epses_value <- epses[5]
#r_value <- thresholds[2]

# Main nested loops 

for (organ in organ_names) {  # Bracket 1 start
  organ_specific_file_list <- file_list[[organ]]  
  for (file_path in organ_specific_file_list) {   # Bracket 2 start
    for (epses_value in epses) { # Bracket 3 start
      for (r_value in thresholds) {               # Bracket 4 start 
        # Track parallel workers
        cat("Using", foreach::getDoParWorkers(), "parallel workers\n")
        
       message("Starting processing combination: ", 
            "Organ:", organ, 
            "File:", basename(file_path), 
            "R2:", r_value, "\n")
        
        match <- str_match(basename(file_path), "CClock_methylation_matrix_?(min\\d+|nofilter)?")
        seqmonk_filter <- match[2]
        print(c(seqmonk_filter))
        
        #------------------------------------------------------
        # STEP 3.1: Load and process methylation data
        #------------------------------------------------------
        
        message("Processing STEP 3.1 - Organ: ", organ, " | Epses: ", epses_value, " | R²: ", r_value)
        
        methylation_data_raw <- data.table::fread(file_path)
        colnames(methylation_data_raw) <- sub("_.*", "", colnames(methylation_data_raw))
        
        methylation_data_raw <- methylation_data_raw %>%
          mutate(Probe = paste(chr, start, end, sep = ":")) %>%
          dplyr::rename("Chromosome"= "chr",
                        "Start" = "start",
                        "End" = "end",
                        "Strand" = "strand")
        
        sample_data_filtered <- sample_data %>%
          dplyr::filter(Sanger_new_ID %in% colnames(methylation_data_raw)) %>%
          as.data.frame()
        
        meth_granges <- GRanges(
          seqnames = methylation_data_raw$Chromosome,
          ranges = IRanges(start = methylation_data_raw$Start, end = methylation_data_raw$End)
        )
        
        dbscan_data <- read.delim(paste0(eps_cluster_dir, "/", organ, "/", seqmonk_filter, "/", r_value, "/CClock_methylation_matrix_", seqmonk_filter, "_clustered_eps_", r_value, "_", epses_value, ".tsv"))
        
        dbscan_granges <- GRanges(
          seqnames = dbscan_data$chr,
          ranges = IRanges(start = dbscan_data$Start, end = dbscan_data$End)
        )
        
        mcols(dbscan_granges)$cluster <- dbscan_data$Cluster
        
        overlaps <- findOverlaps(meth_granges, dbscan_granges)
        
        methylation_data_clustered <- methylation_data_raw[overlaps@from,]
        methylation_data_clustered$cluster <- dbscan_granges$cluster[overlaps@to]
        
        methylation_data_clustered$unique_ID <- paste(methylation_data_clustered$Chromosome, methylation_data_clustered$cluster, sep = "_")
        
        cluster_id <- methylation_data_clustered[, c("Probe","unique_ID")]
        
        id_columns <- c("Chromosome", "Start", "End", "Strand", "Probe", "cluster", "unique_ID")
        value_columns <- setdiff(names(methylation_data_raw), id_columns)

        methylation_data_clustered_mean <- methylation_data_clustered %>%
          group_by(unique_ID) %>%
select(all_of(value_columns)) %>%
          summarise(
            across(
              .fns = ~ mean(.x, na.rm = TRUE)
            )
          )
        
        methylation_data_long <- data.table::melt(
          data.table::as.data.table(methylation_data_clustered),
          id.vars = id_columns,
          measure.vars = value_columns,
          variable.name = "Sanger_new_ID",
          value.name = "Methylation_level"
        ) %>%
          dplyr::left_join(sample_data_filtered, by = "Sanger_new_ID")
        
        #------------------------------------------------------
        # STEP 3.2: Compute R² and Mean Methylation per Probe
        #------------------------------------------------------
        message("Processing STEP 3.2 - Organ: ", organ, " | Epses: ", epses_value, " | R²: ", r_value)
        if (!is.null(methylation_data_long) && ncol(methylation_data_long) > 0) {
          calculated_metrics <- methylation_data_long %>%
            group_by(Probe) %>%
            summarise(
              r_squared = tryCatch({
                summary(lm(Methylation_level ~ Age, data = cur_data()))$r.squared
              }, error = function(e) NA_real_),
              mean_methylation = mean(Methylation_level, na.rm = TRUE),
              .groups = "drop"
            ) %>%
            as.data.frame()
        } else {
          warning("Missing or incomplete data for = ", organ, ",  at r threshold = ", r_value)
          calculated_metrics <- NULL
        }
        
        if (!is.null(calculated_metrics)) {
          if (r_value == "nofilter") {
            filtered_metrics <- calculated_metrics
          } else {
            filtered_metrics <- calculated_metrics %>%
              filter(!is.na(r_squared) & r_squared > as.numeric(r_value))
          }
        } else {
          warning("Missing or incomplete data for = ", organ, ",  at r threshold = ", r_value)
          filtered_metrics <- data.frame()
        }
        
        # Initialize best metrics and model containers
        best_metrics <- list(
          best_train_mae = Inf,
          best_test_mae = Inf,
          best_test_r_squared = -Inf,
          best_train_r_squared = -Inf,
          best_outter_seed = NULL,
          best_inner_seed = NULL,
          best_model = NULL,
          best_nb_cpg = NA,
          best_cpg_ID = NULL,
          best_probe_list = NULL,
          best_val_mae = Inf,
          best_val_r_squared = -Inf
        )
        
        failed_models <- list()
        validation_metrics <- list()
        best_models <- list()
        
        #------------------------------------------------------
        # STEP 4: Elastic Net Regression Nested Loops
        #------------------------------------------------------
        message("Processing STEP 4 - Organ: ", organ, " | Epses: ", epses_value, " | R²: ", r_value)
        if (!is.null(filtered_metrics) && nrow(filtered_metrics) > 1) {
          current_probes <- filtered_metrics
        } else {
          warning("Missing or incomplete data for = ", organ, ",  at r threshold = ", r_value)
          current_probes <- data.frame()
        }
        
        for (i in 1:num_runs) { # Outer loop start
          tryCatch({
            message(paste0("Starting outer loop ", i, " R threshold = ", r_value))
            set.seed(i)
            outerloop_current_seed <- i
            
            # Split data for clock creation and final validation
            if (!is.null(sample_data_filtered) && all(c("Age", "Organ") %in% colnames(sample_data_filtered))) {
              validationIndex <- createDataPartition(
                y = paste(sample_data_filtered$Age, sample_data_filtered$Organ),
                p = 0.8,
                list = FALSE,
                times = 1
              )
            } else {
              warning("Missing or incomplete data for = ", organ, ",  at r threshold = ", r_value)
              validationIndex <- numeric(0)
            }
            
            if (!is.null(validationIndex)) {
              final_validation_metadata <- sample_data_filtered[-validationIndex, ]
              creation_clock_metadata <- sample_data_filtered[validationIndex, ]
            } else {
              warning("Missing or incomplete data for = ", organ, ",  at r threshold = ", r_value)
              final_validation_metadata <- data.frame()
              creation_clock_metadata <- data.frame()
            }
            
            # Inner loop for CpG selection and model training
            for (j in 1:num_runs) {
              tryCatch({
                message(paste0("Starting outer loop: ", i, ", inner loop: ", j, ", R threshold = ", r_value))
                set.seed(j)
                innerloop_current_seed <- j
                
                # Randomly select number and subset CpGs
                if (!is.null(current_probes) && nrow(current_probes) > 1) {
                  num_selected_cpgs <- sample(1:nrow(current_probes), 1)
                  selected_cpgs <- current_probes[sample(nrow(current_probes), num_selected_cpgs), , drop = FALSE]
                  
                  methylation_data_filtered <- methylation_data_clustered %>%
                    filter(Probe %in% selected_cpgs$Probe)
                  
                  methylation_data_filtered <- column_to_rownames(methylation_data_filtered, "Probe")
                  methylation_data_filtered <- methylation_data_filtered %>%
                    dplyr::select(-Chromosome, -Start, -End, -Strand, -cluster, -unique_ID)
                } else {
                  warning("Missing or incomplete data for = ", organ, ",  at r threshold = ", r_value)
                  num_selected_cpgs <- numeric(0)
                  selected_cpgs <- numeric(0)
                  methylation_data_filtered <- data.frame()
                }
                
                # Data cleaning
                if (!is.null(methylation_data_filtered) && ncol(methylation_data_filtered) > 0) {
                  methylation_data_filtered <- methylation_data_filtered[rowSums(is.na(methylation_data_filtered)) != ncol(methylation_data_filtered), ]
                  
                  # Check IDs match
                  print(all(sample_data_filtered$Sanger_new_ID %in% colnames(methylation_data_filtered)))
                  
                  methylation_data_filtered1 <- na.omit(methylation_data_filtered)
                  
                  all_zeros_index <- colSums(methylation_data_filtered, na.rm = T) == 0
                  methylation_data_filtered <- methylation_data_filtered[, !all_zeros_index]
                  
                  all_zero_rows <- rowSums(methylation_data_filtered[, sapply(methylation_data_filtered, is.numeric)], na.rm = TRUE) == 0
                  methylation_data_filtered <- methylation_data_filtered[!all_zero_rows, ]
                } else {
                  warning("Missing or incomplete data for = ", organ, ",  at r threshold = ", r_value)
                  methylation_data_filtered <- data.frame()
                }
                
                # Subset metadata and methylation data for training
                if (!is.null(creation_clock_metadata) && nrow(creation_clock_metadata) > 0 && nrow(methylation_data_filtered) > 0) {
                  rownames(creation_clock_metadata) <- creation_clock_metadata$Sanger_new_ID
                  
                  valid_ids <- intersect(rownames(creation_clock_metadata), colnames(methylation_data_filtered))
                  methylation_data_filtered_final <- methylation_data_filtered[, valid_ids, drop = FALSE]
                  methylation_data_filtered_final <- t(methylation_data_filtered_final)
                  
                  creation_clock_metadata <- creation_clock_metadata[order(rownames(creation_clock_metadata)), ]
                  methylation_data_filtered_final <- methylation_data_filtered_final[order(rownames(methylation_data_filtered_final)), ]
                  
                  if (!all(rownames(creation_clock_metadata) == rownames(methylation_data_filtered_final))) {
                    stop("Error: Row names are not aligned after sorting!")
                  }
                } else {
                  warning("Missing or incomplete data for = ", organ, ",  at r threshold = ", r_value)
                  methylation_data_filtered_final <- data.frame()
                }
                
                # Split training and testing data from creation clock metadata
                if (!is.null(creation_clock_metadata) && ncol(methylation_data_filtered_final) > 1 && is.matrix(methylation_data_filtered_final)) {
                  trainIndex.all <- createDataPartition(
                    y = paste(creation_clock_metadata$Age, creation_clock_metadata$Organ), 
                    p = 0.75, list = FALSE, times = 1
                  )
                  
                  metadata_train.final <- creation_clock_metadata[trainIndex.all, ]
                  metadata_test.final <- creation_clock_metadata[-trainIndex.all, ]
                  
                  methylation_data_train.final <- methylation_data_filtered_final[trainIndex.all, ]
                  methylation_data_test.final <- methylation_data_filtered_final[-trainIndex.all, ]
                  
                  if (!all(rownames(methylation_data_train.final) == rownames(metadata_train.final))) {
                    stop("Error: Row names are not aligned after sorting!")
                  }
                } else {
                  warning("Missing or incomplete data for = ", organ, ",  at r threshold = ", r_value)
                  metadata_train.final <- data.frame()
                  metadata_test.final <- data.frame()
                  methylation_data_train.final <- data.frame()
                  methylation_data_test.final <- data.frame()
                }
                
                # Remove columns with NAs from train and subset test accordingly
                if (!is.null(methylation_data_train.final) && nrow(methylation_data_train.final) > 0) {
                  methylation_data_train.final <- methylation_data_train.final[, colSums(is.na(methylation_data_train.final)) != nrow(methylation_data_train.final)]
                  methylation_data_test.final <- methylation_data_test.final[, colnames(methylation_data_test.final) %in% colnames(methylation_data_train.final)]
                } else {
                  warning("Missing or incomplete data for = ", organ, ",  at r threshold = ", r_value)
                  methylation_data_train.final <- data.frame()
                  methylation_data_test.final <- data.frame()
                }
                
                # Model fitting and prediction
                if (ncol(methylation_data_train.final) > 1) {
                  control <- trainControl(
                    method = "repeatedcv",
                    number = 10,
                    repeats = 50,
                    search = "random",
                    verboseIter = TRUE,
                    allowParallel = TRUE
                  )
                  


                  elastic_model_train <- caret::train(
                    x = methylation_data_train.final,
                    y = metadata_train.final$Age,
                    method = "glmnet",
                    tuneLength = 25,
                    trControl = control,
                    preProcess = "medianImpute" # This replaces all NAs with the median of that feature
                  )
                  
                  coefficients <- broom::tidy(elastic_model_train$finalModel, s = elastic_model_train$bestTune$lambda)
                  selected_coeffs <- coefficients[coefficients$estimate != 0, ]
                  
                  predict_age <- function(modell, data) {
                    predicted_values <- unname(predict(modell, newdata = data))
                    predicted_values[predicted_values < 0] <- 0
                    data.frame(Predicted_Age = predicted_values)
                  }
                  
                  train_predicted_ages <- predict_age(elastic_model_train, methylation_data_train.final)
                  test_predicted_ages <- predict_age(elastic_model_train, methylation_data_test.final)
                  
                  train_mae <- mean(abs(train_predicted_ages$Predicted_Age - metadata_train.final$Age))
                  test_mae <- mean(abs(test_predicted_ages$Predicted_Age - metadata_test.final$Age))
                  train_r_squared <- cor(train_predicted_ages$Predicted_Age, metadata_train.final$Age)^2
                  test_r_squared <- cor(test_predicted_ages$Predicted_Age, metadata_test.final$Age)^2
                  
                } else {
                  warning("Missing or incomplete data for = ", organ, ",  at r threshold = ", r_value)
                  train_mae <- NA_real_
                  test_mae <- NA_real_
                  train_r_squared <- NA_real_
                  test_r_squared <- NA_real_
                }
                
                # Keep only best model by test R squared
                if (!is.na(test_r_squared) && test_r_squared > best_metrics$best_test_r_squared) {
                  best_metrics$best_train_mae <- train_mae
                  best_metrics$best_test_mae <- test_mae
                  best_metrics$best_train_r_squared <- train_r_squared
                  best_metrics$best_test_r_squared <- test_r_squared
                  best_metrics$best_outter_seed <- outerloop_current_seed
                  best_metrics$best_inner_seed <- innerloop_current_seed
                  best_metrics$best_nb_cpg <- num_selected_cpgs
                  best_metrics$best_cpg_ID <- selected_cpgs$Probe
                  best_metrics$best_model <- elastic_model_train
                  best_metrics$best_probe_list <- methylation_data_filtered_final
                }
                
              }, error = function(e) {
                message("Inner loop error: ", e$message)
                failed_models[[length(failed_models) + 1]] <<- list(
                  outerloop = i,
                  innerloop = j,
                  organ = organ,
                  r_value = r_value,
                  error = e$message
                )
              })
            } # Inner loop end
            
            message(paste0("Finished inner loop ", j))
          }, error = function(e) {
            message("Outer loop error: ", e$message)
          })
        } # Outer loop end
        
        # Save best model and results
        organ_output_path <- file.path(output_directory, organ, r_value)
        dir.create(organ_output_path, recursive = TRUE, showWarnings = FALSE)
        
        if (!is.null(best_metrics$best_model)) {
          saveRDS(best_metrics$best_model, file = file.path(organ_output_path, paste0("best_model_", organ, "_r_", r_value, ".rds")))
        }
        
        # Save metrics
        metrics_file <- file.path(organ_output_path, paste0("best_metrics_", organ, "_r_", r_value, ".rds"))
        saveRDS(best_metrics, metrics_file)
        
        message(paste0("Completed R threshold: ", r_value, " for organ: ", organ))
        
      } # Bracket 4 end (threshold loop)
    }   # Bracket 3 end (eps loop)
  }     # Bracket 2 end (file loop)
}       # Bracket 1 end (organ loop)

stopCluster(cl)
message("All training complete.")
