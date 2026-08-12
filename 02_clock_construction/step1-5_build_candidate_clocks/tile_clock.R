# ==========================================================================
# PIPELINE STAGE: 02 / Clock construction - steps 1-5 (Tile)
# ==========================================================================
# Builds candidate Tile epigenetic clocks: data prep, probe filtering,
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



# Tile Epigenetic Clock Construction Pipeline (Steps 1–6)
#
# This pipeline processes DNA methylation data across multiple organs to build point epigenetic clocks, that is based on tiled GPG sites form both whole genome or published . 
# Step 1–2: Data preparation and probe filtering based on quality and coverage.
# Step 3: Computation of per-probe R² and mean methylation metrics, applying R² thresholds.
# Step 4: Initialise modelling scaffolds per organ and R² threshold using filtered CpGs.
# Step 5: Filter CpG probe sets by R² threshold and prepare data structures for modelling.
# Step 6: Iterative lasso regression model training, testing, and independent validation 
#         with random CpG subset selection and best model selection based on prediction accuracy.
#
# To optimise memory usage, outputs from each step are saved to disk, allowing incremental processing 
# without requiring all data to be held in memory simultaneously.
#
# The pipeline outputs model objects, filtered probe lists, and performance metrics organised 
# by organ and R² threshold, facilitating robust age prediction model construction.

# ======================
# Setup and Dependencies
# ======================
required_packages <- c(
  "data.table", "dplyr", "tidyr", "readxl", "stringr", "dbscan",
  "progress", "caret", "glmnetUtils", "doParallel", "shades",
  "RColorBrewer", "glmnet", "tidyverse", "gridExtra", "patchwork",
  "progressr", "purrr", "rlang"
)

install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cran.rstudio.com/")
  }
}
invisible(lapply(required_packages, install_if_missing))
invisible(lapply(required_packages, function(pkg) {
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}))

# ======================
# Parallel Setup
# ======================
n_cores <- parallel::detectCores() - 1
cl <- makeCluster(n_cores)
registerDoParallel(cl)


# ======================
# Define parameters and directories
# ======================
input_directory_base <- file.path(PROJECT_ROOT, "multitissue/making_annotation_files/R/Results/tile_clock")
input_directory_metadata <- file.path(PROJECT_ROOT, "multitissue")
output_directory <- file.path(PROJECT_ROOT, "multitissue/making_clocks/R/results/tile_clock_final_final")

dir.create(output_directory, showWarnings = FALSE, recursive = TRUE)

# ======================
# Define R-squared thresholds
# ======================
thresholds <- c(
  "r2_nofilter" = "nofilter",
  "r2_0.05" = "0.05",
  "r2_0.1" = "0.1",
  "r2_0.2" = "0.2",
  "r2_0.3" = "0.3"
)

# ======================
# Load metadata
# ======================
cat("Loading metadata...\n")
sample_data <- read_excel(file.path(input_directory_metadata, "final_metadata_uniqueID_new.xlsx")) %>%
  dplyr::select(Clock, sequencing_block, plate_name, Sanger_new_ID, Organ, original_tube_ID, Age, Sex, Treatment) %>%
  filter(Clock == "multitissue")

cat("Metadata loaded: ", nrow(sample_data), " rows\n")

# ======================
# Define ID columns
# ======================
id_columns <- c("Probe", "chr", "start", "end","strand")

# ======================
# Function to process methylation data files in a folder
# ======================
process_organ_folder_long <- function(organ_folder_path, organ_name, tile_type_name, r2_name, output_base_dir) {
  message("Processing STEP 2 - Organ: ", organ_name, " | Tile: ", tile_type_name, " | R²: ", r2_name)
  
  file.list <- list.files(path = organ_folder_path, pattern = "\\.tsv$", full.names = TRUE)
  if (length(file.list) == 0) {
    warning("No files found in ", organ_folder_path)
    return(NULL)
  }
  
  file.names <- tools::file_path_sans_ext(basename(file.list))
  
  for (i in seq_along(file.list)) {
    f <- file.list[[i]]
    fname <- file.names[[i]]
    
    df <- fread(f)
    colnames(df) <- sub("_.*", "", colnames(df))
    df <- df %>% mutate(Probe = paste(chr, start, end, sep = ":"))
    
    sample_data_filtered <- sample_data %>%
      filter(Sanger_new_ID %in% colnames(df)) %>%
      as.data.frame()
    
    df_long <- melt(
      as.data.table(df),
      id.vars = id_columns,
      variable.name = "Sanger_new_ID",
      value.name = "Methylation_level"
    ) %>%
      left_join(sample_data_filtered, by = "Sanger_new_ID")
    
    # Corrected path: tile -> organ -> r²
    save_path <- file.path(output_base_dir, "step2", tile_type_name, organ_name, r2_name)
    dir.create(save_path, recursive = TRUE, showWarnings = FALSE)
    
    # Save .rds
    saveRDS(df_long, file = file.path(save_path, paste0(fname, "_long.rds")))
    
    # Clean up
    rm(df, df_long, sample_data_filtered)
    gc()
  }
}

# ================
# STEP 2 Loop: All Tiles, Organs, and R² thresholds
# ================
cat("Starting main processing loops...\n")

# List tile types (one level below input_directory_base)
tile_types <- list.dirs(path = input_directory_base, recursive = FALSE, full.names = TRUE)
tile_type_names <- basename(tile_types)

cat("Tile types found: ", paste(tile_type_names, collapse = ", "), "\n")

# Cache organ directories and names for each tile type
organ_dirs_cache <- vector("list", length(tile_types))
organ_names_cache <- vector("list", length(tile_types))

for (j in seq_along(tile_types)) {
  organ_dirs_cache[[j]] <- list.dirs(path = tile_types[j], recursive = FALSE, full.names = TRUE)
  organ_names_cache[[j]] <- basename(organ_dirs_cache[[j]])
  cat("Tile type '", tile_type_names[j], "' has organs: ", paste(organ_names_cache[[j]], collapse = ", "), "\n")
}

for (r2_name in names(thresholds)) {
  cat("Processing R-squared threshold: ", r2_name, "\n")
  
  for (j in seq_along(tile_types)) {
    tile_type_path <- tile_types[j]
    tile_type_name <- tile_type_names[j]
    
    organ_dirs <- organ_dirs_cache[[j]]
    organ_names <- organ_names_cache[[j]]
    
    for (i in seq_along(organ_dirs)) {
      organ_dir <- organ_dirs[i]
      organ_name <- organ_names[i]
      
      cat("Processing organ: ", organ_name, " (Tile: ", tile_type_name, ")\n")
      
      # Run step 2 for this combination
      try({
        process_organ_folder_long(
          organ_folder_path = organ_dir,
          organ_name = organ_name,
          tile_type_name = tile_type_name,
          r2_name = r2_name,
          output_base_dir = output_directory
        )
      }, silent = TRUE)
    }
  }
}

message("Step 2 complete — files saved in: step2/<tile>/<organ>/<r2>/")


# ================
# STEP 3 - Compute R² and Mean Methylation from saved .rds files
# ================
# Directories for input/output
input_root <- file.path(output_directory, "step2")
output_root <- file.path(output_directory, "step3")
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

# Iterate over thresholds and tile types
for (r2_name in names(thresholds)) {
  r2_threshold <- thresholds[[r2_name]]
  
  for (tile_type_name in tile_type_names) {
    tile_dir <- file.path(input_root, tile_type_name)
    if (!dir.exists(tile_dir)) next
    
    organ_dirs <- list.dirs(tile_dir, recursive = FALSE)
    
    for (organ_dir in organ_dirs) {
      organ_name <- basename(organ_dir)
      organ_threshold_dir <- file.path(organ_dir, r2_name)
      if (!dir.exists(organ_threshold_dir)) next
      
      rds_files <- list.files(organ_threshold_dir, pattern = "_long\\.rds$", full.names = TRUE)
      if (length(rds_files) == 0) next
      
      for (rds_file in rds_files) {
        sample_name <- sub("_long\\.rds$", "", basename(rds_file))
        message("Processing STEP 3 : ", tile_type_name, " | ", organ_name, " | ", r2_name, " | ", sample_name)
        
        df <- readRDS(rds_file)
        if (!all(c("Probe", "Methylation_level", "Age") %in% colnames(df))) {
          warning("Required columns missing in: ", rds_file)
          next
        }
        
        metrics_df <- df %>%
          group_by(Probe) %>%
          summarise(
            r_squared = tryCatch({
              summary(lm(Methylation_level ~ Age, data = cur_data()))$r.squared
            }, error = function(e) NA_real_),
            mean_methylation = mean(Methylation_level, na.rm = TRUE),
            .groups = "drop"
          )
        
        # Construct output directory path
        metrics_out_dir <- file.path(output_root, tile_type_name, organ_name, r2_name)
        dir.create(metrics_out_dir, recursive = TRUE, showWarnings = FALSE)
        
        # Save metrics
        saveRDS(metrics_df, file = file.path(metrics_out_dir, paste0(sample_name, "_metrics.rds")))
        
        # Optional filtering step
        if (r2_threshold != "nofilter") {
          filtered_df <- metrics_df %>%
            filter(!is.na(r_squared) & r_squared > as.numeric(r2_threshold))
          
          saveRDS(filtered_df, file = file.path(metrics_out_dir, paste0(sample_name, "_metrics.rds")))
        }
        
        # Memory cleanup
        rm(df, metrics_df, filtered_df)
        gc()
      }
    }
  }
}


message("Step 3 complete — metrics and filtered metrics saved by organ and threshold.")


#------------------------------------------------------
# STEP 4: Setup Epigenetic Clock Model Infrastructure
#------------------------------------------------------

# Define the age prediction function
predict_age <- function(x, y) {
  predicted_values <- predict(x, newx = y)
  colnames(predicted_values) <- "Predicted_Age"
  return(predicted_values)
}

# Input and output directories
input_dir <- file.path(output_directory, "step3")
output_dir <- file.path(output_directory, "step4")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# List all tile types in step3
tile_dirs <- list.dirs(input_dir, recursive = FALSE, full.names = TRUE)

for (tile_dir in tile_dirs) {
  tile_type <- basename(tile_dir)
  
  # List organs for this tile
  organ_dirs <- list.dirs(tile_dir, recursive = FALSE, full.names = TRUE)
  
  for (organ_dir in organ_dirs) {
    organ_name <- basename(organ_dir)
    
    # List R² directories for this organ
    r2_dirs <- list.dirs(organ_dir, recursive = FALSE, full.names = TRUE)
    
    for (r2_dir in r2_dirs) {
      r_squared_value <- basename(r2_dir)
      
      # Define output directory: tile → organ → r2
      output_r2_dir <- file.path(output_dir, tile_type, organ_name, r_squared_value)
      dir.create(output_r2_dir, recursive = TRUE, showWarnings = FALSE)
      
      # List all filtered metrics files
      metric_files <- list.files(r2_dir, pattern = "_metrics\\.rds$", full.names = TRUE)
      
      for (file in metric_files) {
        file_name <- basename(file)
        sample_id <- str_remove(file_name, "_metrics\\.rds$")
        
        message("Processing STEP 4: ", tile_type, " | ", organ_name, " | ", r_squared_value, " | ", sample_id)
        
        df_filtered <- readRDS(file)
        
        # Initialise model scaffold
        file_entry <- list(
          filtered_metrics = df_filtered,
          best_metrics = list(
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
          ),
          failed_models = list(),
          validation_metrics = list(),
          best_models = list()
        )
        
        # Save model scaffold
        output_file <- file.path(output_r2_dir, paste0(sample_id, "_clock_setup.rds"))
        saveRDS(file_entry, output_file, compress = "xz")
      }
    }
  }
}

message("✅ Step 4 completed: Initialised modelling scaffolds under step4/<tile>/<organ>/<r2>/")


#------------------------------------------------------
# STEP 5: Extract filtered CpGs from Step 4 saved files per organ and R² threshold
#------------------------------------------------------

# Input/output directories
input_dir <- file.path(output_directory, "step4")
output_dir <- file.path(output_directory, "step5")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

message("Starting STEP 5 - extracting filtered CpGs per tile type, organ and R² threshold")

# List all tile-type directories (e.g., '50bp', '1kb', etc.)
tile_dirs <- list.dirs(input_dir, recursive = FALSE, full.names = TRUE)

for (tile_dir in tile_dirs) {
  tile_type <- basename(tile_dir)
  
  # List organ folders under each tile type
  organ_dirs <- list.dirs(tile_dir, recursive = FALSE, full.names = TRUE)
  
  for (organ_dir in organ_dirs) {
    organ_name <- basename(organ_dir)
    
    # List R² threshold folders
    r2_dirs <- list.dirs(organ_dir, recursive = FALSE, full.names = TRUE)
    
    for (r2_dir in r2_dirs) {
      r_squared_value <- basename(r2_dir)
      
      # Define output directory: tile → organ → r2
      output_r2_dir <- file.path(output_dir, tile_type, organ_name, r_squared_value)
      dir.create(output_r2_dir, recursive = TRUE, showWarnings = FALSE)
      
      # List all clock setup files (from Step 4)
      clock_setup_files <- list.files(r2_dir, pattern = "_clock_setup\\.rds$", full.names = TRUE)
      
      for (file in clock_setup_files) {
        file_name <- basename(file)
        sample_id <- str_remove(file_name, "_clock_setup\\.rds$")
        
        message("Processing STEP 5: ", tile_type, " | ", organ_name, " | ", r_squared_value, " | ", sample_id)
        
        file_entry <- readRDS(file)
        
        # Extract filtered_metrics and assign to current_probes if valid
        if (!is.null(file_entry$filtered_metrics) && nrow(file_entry$filtered_metrics) > 1) {
          current_probes <- file_entry$filtered_metrics
        } else {
          warning("Missing or incomplete filtered metrics for ", sample_id,
                  " at R² threshold ", r_squared_value)
          current_probes <- data.frame()
        }
        
        # Save per sample
        output_file <- file.path(output_r2_dir, paste0(sample_id, "_current_probes.rds"))
        saveRDS(current_probes, output_file)
      }
    }
  }
}

message("✅ Step 5 completed: filtered CpG probe lists saved under step5/<tile>/<organ>/<r2>/")


#------------------------------------------------------
# Step 6: Epigenetic Clock Model Training, Validation, and Selection
#         Across Tiles, Organs, and R² Thresholds
#------------------------------------------------------

# Parameters
num_runs <- 1000  # 
input_dir <- file.path(output_directory, "step4") # input: clock scaffoldsform step 4 
input_dir_2 <- file.path(output_directory, "step2") # input: clock large metadata step 2
output_dir <- file.path(output_directory, "step6") # output: updated scaffolds with models
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)


# Define the age prediction function
predict_age <- function(x, y) {
  predicted_values <- predict(x, newx = y)
  colnames(predicted_values) <- "Predicted_Age"
  return(predicted_values)
}

# Get tile folders
tile_dirs <- list.dirs(input_dir, recursive = FALSE, full.names = TRUE)

for (tile_dir in tile_dirs) {
  tile_name <- basename(tile_dir)
  
  # Get organ folders inside tile
  organ_dirs <- list.dirs(tile_dir, recursive = FALSE, full.names = TRUE)
  
  for (organ_dir in organ_dirs) {
    organ_name <- basename(organ_dir)
    
    # Get R² threshold folders inside organ
    r2_dirs <- list.dirs(organ_dir, recursive = FALSE, full.names = TRUE)
    
    for (r2_dir in r2_dirs) {
      r_squared_value <- basename(r2_dir)
      
      # Output directory for this tile-organ-R2 combination
      output_r2_dir <- file.path(output_dir, tile_name, organ_name, r_squared_value)
      dir.create(output_r2_dir, recursive = TRUE, showWarnings = FALSE)
      
      # List sample scaffold files
      scaffold_files <- list.files(r2_dir, pattern = "_clock_setup\\.rds$", full.names = TRUE)
      
      for (scaffold_file in scaffold_files) {
        sample_id <- str_remove(basename(scaffold_file), "_clock_setup\\.rds$")
        
        message("Processing STEP 6: ", tile_name, " | ", organ_name, " | ", r_squared_value, " | ", sample_id)
        
        file_entry <- readRDS(scaffold_file)
        file_entry$methylation_data <- fread(paste0(input_directory_base,"/",tile_name,"/",organ_name,"/",sample_id,".tsv"))
        colnames(file_entry$methylation_data) <- sub("_.*", "", colnames(file_entry$methylation_data))
        file_entry$methylation_data <- file_entry$methylation_data %>% mutate(Probe = paste(chr, start, end, sep = ":"))
        rownames(file_entry$methylation_data) <- file_entry$methylation_data$Probe
        file_entry$sample_data_filtered <- sample_data %>%
          filter(Sanger_new_ID %in% colnames(file_entry$methylation_data)) %>%
          as.data.frame()
        
        id_columns <- c("Probe", "chr", "start", "end","strand")
        # Safely wrap the entire modelling & validation in tryCatch
        tryCatch({
          for (i in 1:num_runs) {
            message(paste0("Outer loop run ", i, " for R2 threshold ", r_squared_value))
            set.seed(i)
            outerloop_current_seed <- i
            
            # 6.1: Create training/validation split (80/20)
            if (!is.null(file_entry$sample_data_filtered) &&
                all(c("Age", "Organ") %in% colnames(file_entry$sample_data_filtered))) {
              
              validationIndex <- createDataPartition(
                y = paste(file_entry$sample_data_filtered$Age, file_entry$sample_data_filtered$Organ),
                p = 0.8,
                list = FALSE,
                times = 1
              )
              file_entry$validationIndex <- validationIndex
              
              # Final validation metadata (20%)
              file_entry$final_validation_metadata <- file_entry$sample_data_filtered[-validationIndex, ]
              # Creation clock metadata (80%)
              file_entry$creation_clock_metadata <- file_entry$sample_data_filtered[validationIndex, ]
            } else {
              warning("Missing or incomplete sample_data_filtered for sample ", sample_id)
              next
            }
            
            for (j in 1:num_runs) {
              message(paste0("Inner loop run ", j, " for R2 threshold ", r_squared_value))
              set.seed(j)
              innerloop_current_seed <- j
              
              if (!is.null(file_entry$current_probes) && nrow(file_entry$current_probes) > 1) {
                num_selected_cpgs <- sample(1:nrow(file_entry$current_probes), 1)
                selected_cpgs <- file_entry$current_probes[sample(nrow(file_entry$current_probes), num_selected_cpgs), , drop = FALSE]
                methylation_data_filtered <- file_entry$methylation_data[rownames(file_entry$methylation_data) %in% selected_cpgs$Probe, , drop = FALSE]
                rownames(methylation_data_filtered) <- methylation_data_filtered$Probe
                
                file_entry$num_selected_cpgs <- num_selected_cpgs
                file_entry$selected_cpgs <- selected_cpgs
                file_entry$methylation_data_filtered <- methylation_data_filtered
              } else {
                warning("Missing or incomplete current_probes for sample ", sample_id)
                next
              }
              
              # 6.3 Clean methylation data
              if (!is.null(file_entry$methylation_data_filtered) &&
                  ncol(file_entry$methylation_data_filtered) > 0) {
                md_filtered <- file_entry$methylation_data_filtered
                md_filtered <- md_filtered[rowSums(is.na(md_filtered)) != ncol(md_filtered), ]
                md_filtered1 <- na.omit(md_filtered)
                rownames(md_filtered) <- md_filtered$Probe
                md_filtered <- md_filtered %>% select(-all_of(id_columns))
                file_entry$methylation_data_filtered <- md_filtered
              } else {
                warning("Incomplete methylation_data_filtered for sample ", sample_id)
                next
              }
              
              # 6.4 Create final filtered matrix & metadata
              if (!is.null(file_entry$creation_clock_metadata) &&
                  nrow(file_entry$creation_clock_metadata) > 0 &&
                  nrow(file_entry$methylation_data_filtered) > 0) {
                
                rownames(file_entry$creation_clock_metadata) <- file_entry$creation_clock_metadata$Sanger_new_ID
                
                md <- file_entry$methylation_data_filtered
                
                valid_ids <- intersect(rownames(file_entry$creation_clock_metadata), colnames(file_entry$methylation_data_filtered))
                md_final <- file_entry$methylation_data_filtered[, ..valid_ids, drop = FALSE]
                rownames(md_final) <- rownames(md)
                
                # Transpose
                md_final_t <- t(md_final)
                
                # assign  row/col names after transpose
                #  rownames(md_final_t) <- colnames(md_final)
                colnames(md_final_t) <- rownames(md_final) 
                
                
                
                file_entry$creation_clock_metadata <- file_entry$creation_clock_metadata[order(rownames(file_entry$creation_clock_metadata)), ]
                md_final_t <-  md_final_t[order(rownames( md_final_t)), ]
                
                if (!all(rownames(file_entry$creation_clock_metadata) == rownames(md_final_t))) {
                  stop("Row names misalignment after sorting for sample ", sample_id)
                }
                
                file_entry$methylation_data_filtered_final <-  md_final_t
              } else {
                warning("Missing metadata or methylation data for sample ", sample_id)
                next
              }
              
              # 6.5 Split creation_clock_metadata into train/test 75%/25%
              if (!is.null(file_entry$creation_clock_metadata) &&
                  ncol(file_entry$methylation_data_filtered_final) > 1 &&
                  is.matrix(file_entry$methylation_data_filtered_final)) {
                
                trainIndex.all <- createDataPartition(
                  y = paste(file_entry$creation_clock_metadata$Age, file_entry$creation_clock_metadata$Organ),
                  p = 0.75,
                  list = FALSE,
                  times = 1
                )
                
                file_entry$metadata_train.final <- file_entry$creation_clock_metadata[trainIndex.all, ]
                file_entry$metadata_test.final <- file_entry$creation_clock_metadata[-trainIndex.all, ]
                file_entry$methylation_data_train.final <- file_entry$methylation_data_filtered_final[trainIndex.all, , drop = FALSE]
                file_entry$methylation_data_test.final <- file_entry$methylation_data_filtered_final[-trainIndex.all, , drop = FALSE]
                
                if (!all(rownames(file_entry$methylation_data_train.final) == rownames(file_entry$metadata_train.final))) {
                  stop("Row names misalignment between train methylation and metadata for sample ", sample_id)
                }
              } else {
                warning("Insufficient data for train/test split for sample ", sample_id)
                next
              }
              
              # 6.6 Remove columns with NA in train and match test columns
              if (nrow(file_entry$methylation_data_train.final) > 0) {
                cols_to_keep <- colSums(is.na(file_entry$methylation_data_train.final)) != nrow(file_entry$methylation_data_train.final)
                file_entry$methylation_data_train.final <- file_entry$methylation_data_train.final[, cols_to_keep, drop = FALSE]
                file_entry$methylation_data_test.final <- file_entry$methylation_data_test.final[, colnames(file_entry$methylation_data_test.final) %in% colnames(file_entry$methylation_data_train.final), drop = FALSE]
              } else {
                warning("Empty train data after NA filtering for sample ", sample_id)
                next
              }
              
              
              # 6.7 Fit lasso model and predict
              if (ncol(file_entry$methylation_data_train.final) > 1) {
                message(paste0("Fitting lasso model for sample ", sample_id))
                
                
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
                
                # Replace the original
                file_entry$methylation_data_train.final <- x_train
                
                # -- for test
                x_test <- file_entry$methylation_data_test.final
                x_test <- as.matrix(x_test)  # ensure it’s a matrix
                
                # Median impute NAs column-wise
                for (i in seq_len(ncol(x_test))) {
                  if (anyNA(x_test[, i])) {
                    x_test[, i][is.na(x_test[, i])] <- median(x_test[, i], na.rm = TRUE)
                  }
                }
                
                # Replace the original
                file_entry$methylation_data_test.final <- x_test
                
                # now making clock
                cv_fit <- cv.glmnet(
                  x = file_entry$methylation_data_train.final,
                  y = file_entry$metadata_train.final$Age,
                  alpha = 1,
                  nfolds = 10,
                  parallel = TRUE,
                  nlambda = 100
                )
                
                lasso_fit <- glmnet(
                  x = file_entry$methylation_data_train.final,
                  y = file_entry$metadata_train.final$Age,
                  lambda = cv_fit$lambda.min
                )
                
                coefficients <- coef(lasso_fit)
                selected_coeffs <- coefficients[, 1] != 0
                
                train_predicted_ages <- predict_age(y = file_entry$methylation_data_train.final, x = lasso_fit)
                test_predicted_ages <- predict_age(y = file_entry$methylation_data_test.final, x = lasso_fit)
                
                train_mae <- mean(abs(train_predicted_ages - file_entry$metadata_train.final$Age))
                test_mae <- mean(abs(test_predicted_ages - file_entry$metadata_test.final$Age))
                
                train_r_squared <- cor(train_predicted_ages, file_entry$metadata_train.final$Age)^2
                test_r_squared <- cor(test_predicted_ages, file_entry$metadata_test.final$Age)^2
                
                file_entry$cv_fit <- cv_fit
                file_entry$lasso_fit <- lasso_fit
                file_entry$coefficients <- coefficients
                file_entry$selected_coeffs <- selected_coeffs
                file_entry$train_predicted_ages <- train_predicted_ages
                file_entry$test_predicted_ages <- test_predicted_ages
                file_entry$train_mae <- train_mae
                file_entry$test_mae <- test_mae
                file_entry$train_r_squared <- train_r_squared
                file_entry$test_r_squared <- test_r_squared
                
              } else {
                warning("Insufficient CpGs for lasso model fitting for sample ", sample_id)
                next
              }
              
              
              # 6.8 Keep only best model by test R²
              if (!is.na(file_entry$test_r_squared) &&
                  file_entry$test_r_squared > file_entry$best_metrics$best_test_r_squared) {
                
                file_entry$best_metrics$best_train_mae <- file_entry$train_mae
                file_entry$best_metrics$best_test_mae <- file_entry$test_mae
                file_entry$best_metrics$best_train_r_squared <- file_entry$train_r_squared
                file_entry$best_metrics$best_test_r_squared <- file_entry$test_r_squared
                file_entry$best_metrics$best_outter_seed <- outerloop_current_seed
                file_entry$best_metrics$best_inner_seed <- innerloop_current_seed
                file_entry$best_metrics$best_nb_cpg <- file_entry$num_selected_cpgs
                file_entry$best_metrics$best_cpg_ID <- file_entry$selected_cpgs
                file_entry$best_metrics$best_probe_list <- r_squared_value
                
                file_entry$best_model <- list(
                  train_mae = file_entry$train_mae,
                  test_mae = file_entry$test_mae,
                  train_r_squared = file_entry$train_r_squared,
                  test_r_squared = file_entry$test_r_squared,
                  cv_fit = file_entry$cv_fit,
                  lasso_fit = file_entry$lasso_fit,
                  metadata_train = file_entry$metadata_train.final,
                  metadata_test = file_entry$metadata_test.final,
                  methylation_data_train = file_entry$methylation_data_train.final,
                  methylation_data_test = file_entry$methylation_data_test.final,
                  train_predicted_ages = file_entry$train_predicted_ages,
                  test_predicted_ages = file_entry$test_predicted_ages
                )
              } else {
                # Store failed model metrics for tracking
                run_name <- paste0("run_outer_", i, "_inner_", j)
                file_entry$failed_models[[run_name]] <- list(
                  nb_cpg = file_entry$num_selected_cpgs,
                  cpg_ID = file_entry$selected_cpgs,
                  train_mae = file_entry$train_mae,
                  test_mae = file_entry$test_mae,
                  train_r_squared = file_entry$train_r_squared,
                  test_r_squared = file_entry$test_r_squared,
                  probe_ID = r_squared_value
                )
              }
            } # end inner loop j
            
            # 6.9 Validation on the original 20% held out data
            if (!is.null(file_entry$best_model) && !is.null(file_entry$final_validation_metadata)) {
              # Filter validation methylation data
              valid_cols <- intersect(
                colnames(file_entry$methylation_data),
                file_entry$final_validation_metadata$Sanger_new_ID
              )
              
              
              
              # Subset only valid columns (plus Probe)
              final_validation_methylation_data <- file_entry$methylation_data[, c("Probe", ..valid_cols), with = FALSE]
              
              
              #  restore probe as row names in a data frame
              final_validation_methylation_data <- as.data.frame(final_validation_methylation_data)
              rownames(final_validation_methylation_data) <- final_validation_methylation_data$Probe
              final_validation_methylation_data$Probe <- NULL
              
              final_validation_methylation_data <- t(final_validation_methylation_data[
                rownames(final_validation_methylation_data) %in% colnames(file_entry$best_model$methylation_data_train),
                , drop = FALSE
              ])
              
              if (ncol(final_validation_methylation_data) > 0) {
                final_val_predicted_ages <- predict_age(y = final_validation_methylation_data, x = file_entry$best_model$lasso_fit)
                
                final_validation_mae <- mean(abs(final_val_predicted_ages - file_entry$final_validation_metadata$Age))
                final_validation_r_squared <- cor(final_val_predicted_ages, file_entry$final_validation_metadata$Age)^2
                
                file_entry$final_validation_model <- list(
                  final_validation_metadata = file_entry$final_validation_metadata,
                  final_validation_methylation_data = final_validation_methylation_data,
                  final_val_predicted_ages = final_val_predicted_ages,
                  final_validation_mae = final_validation_mae,
                  final_validation_r_squared = final_validation_r_squared
                )
              } else {
                warning("No CpGs in validation methylation data for sample ", sample_id)
              }
            }
            
            # 6.10 Keep best validation model
            if (!is.null(file_entry$best_model) && !is.null(file_entry$final_validation_model)) {
              if (!is.na(file_entry$final_validation_model$final_validation_r_squared) &&
                  file_entry$final_validation_model$final_validation_r_squared > file_entry$best_metrics$best_val_r_squared) {
                
                file_entry$best_metrics$best_val_r_squared <- file_entry$final_validation_model$final_validation_r_squared
                file_entry$best_metrics$best_val_mae <- file_entry$final_validation_model$final_validation_mae
                file_entry$best_models$best_validation_metrics <- file_entry$final_validation_model
              }
            }
          } # end outer loop i
          
          # Save updated scaffold with model info
          saveRDS(file_entry, file = file.path(output_r2_dir, paste0(sample_id, "_clock_modelled.rds")), compress = "xz")
          message("Saved modelled scaffold for sample ", sample_id)
        }, error = function(e) {
          message("Error processing sample ", sample_id, ": ", e$message)
        }) # end tryCatch
      } # end sample loop
    } # end R² loop
  } # end organ loop
} # end tile loop

message("Step 6 completed: models fitted and validated for all tiles, organs, thresholds, and samples.")
