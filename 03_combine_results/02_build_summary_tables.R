# ==========================================================================
# PIPELINE STAGE: 03 / Combine results (2 of 2)
# ==========================================================================
# Reads the collated CSVs from 01_collate_model_metrics.R (plus the raw
# .rds model outputs again, for CpG identities) and builds the final
# summary tables and exploratory scatterplots (predicted vs chronological
# age per tissue/model) used as the shared input data for both aim-specific
# downstream analysis notebooks (05_clock_design_benchmarking/,
# 06_cross_tissue_ageing_biology/).
# ==========================================================================

# Set the root directory of this project (edit here, or export EPICLOCK_ROOT
# before launching R / sbatch). All absolute paths below are built from this.
PROJECT_ROOT <- Sys.getenv("EPICLOCK_ROOT", unset = normalizePath("."))

# ---------------------------------------------------------------------------
# Script: Plotting MAE and R-squared from epigenetic clocks & collecting model metrics
# Author: Jean-Charles de Coriolis
#
# Description:
#   Processes epigenetic clock .rds model outputs:
#     1. Extract predicted vs chronological ages (Train, Validation)
#     2. Collect MAE / R² and model metadata
#     3. Generate ggplot2 scatterplots with annotation
#     4. Save combined plots via ggsave
#     5. Incrementally write metrics, predicted ages, CpG sets
# ---------------------------------------------------------------------------

# ------------------- Load Required Packages -------------------
library(dplyr)
library(ggplot2)
library(patchwork)
library(readxl)

# ------------------- Input and Output Paths -------------------
input_folder_C  <- file.path(PROJECT_ROOT, "multitissue/making_clocks/R/results/cluster_clock_new/step4.array.new.new.new")
input_folder_T  <- file.path(PROJECT_ROOT, "multitissue/making_clocks/R/results/tile_clock_final_final/step6.array.new.new")
input_folder_P  <- file.path(PROJECT_ROOT, "multitissue/making_clocks/R/results_split_organ/point_clock/step6.array.new.new")

input_directory_metadata  <- file.path(PROJECT_ROOT, "multitissue")

output_folder_files <- file.path(PROJECT_ROOT, "multitissue/making_clocks/R/scripts/Graphs/results/final_merged")
output_folder_plots <- file.path(PROJECT_ROOT, "multitissue/making_clocks/R/scripts/Graphs/results/final_plots")

dir.create(output_folder_files, recursive = TRUE, showWarnings = FALSE)
dir.create(output_folder_plots, recursive = TRUE, showWarnings = FALSE)

# ------------------- Load Metadata -------------------
metadata_file <- file.path(input_directory_metadata, "final_metadata_uniqueID_new.xlsx")
metadata <- read_excel(metadata_file) %>%
  select(Clock, sequencing_block, plate_name, Sanger_new_ID, Organ,
         original_tube_ID, Age, Sex, Treatment, individual_fish_ID) %>%
  filter(Clock == "multitissue") %>%
  as.data.frame(stringsAsFactors = FALSE)

# ------------------- Palette -------------------
okabe_ito <- c(
  "#E69F00", "#56B4E9", "#009E73", "#F0E442",
  "#0072B2", "#D55E00", "#CC79A7", "#999999"
)
okabe_ito_palette <- colorRampPalette(okabe_ito)

# ------------------- ggplot2 Plotting Function -------------------
plot_predicted_vs_chron_gg <- function(model_list, dataset_name = "Train") {
  
  df <- model_list$data
  df$predicted_age <- as.numeric(model_list$predicted_ages)
  df$Organ <- factor(df$Organ)
  
  mae_value <- round(model_list$mae, 2)
  r2_value  <- round(model_list$r_squared, 2)
  nb_cpg    <- model_list$nb_cpg
  
  organs <- unique(df$Organ)
  cols <- okabe_ito_palette(length(organs))
  names(cols) <- organs
  
  # annotation labels
  lab_mae <- sprintf("MAE: %s", mae_value)
  lab_r2  <- sprintf("R²: %s", r2_value)
  lab_cpg <- if (!is.null(nb_cpg) && dataset_name == "Train") sprintf("CpGs: %s", nb_cpg) else NULL
  
  ymax <- max(df$predicted_age)
  xmin <- min(df$Age)
  
  p <- ggplot(df, aes(x = Age, y = predicted_age, colour = Organ)) +
    geom_point(size = 2, alpha = 0.8) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    scale_colour_manual(values = cols) +
    theme_bw(base_size = 12) +
    labs(
      title = dataset_name,
      x = "Chronological age (months)",
      y = "Predicted age (months)",
      colour = "Organ"
    ) +
    annotate("text", x = xmin, y = ymax, label = lab_mae, hjust = 0, vjust = 1.5, size = 3) +
    annotate("text", x = xmin, y = ymax, label = lab_r2,  hjust = 0, vjust = 3, size = 3)
  
  if (!is.null(lab_cpg)) {
    p <- p +
      annotate("text", x = xmin, y = ymax, label = lab_cpg, hjust = 0, vjust = 4.5, size = 3)
  }
  
  return(p)
}

# ------------------- Safe Write/Append -------------------
append_or_write <- function(df, filepath) {
  if (!file.exists(filepath)) {
    write.table(df, file = filepath, sep = ",", row.names = FALSE, col.names = TRUE)
  } else {
    write.table(df, file = filepath, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
  }
}

# ------------------- Input List -------------------
input_list <- list(
  "P" = input_folder_P,
  "T" = input_folder_T,
  "C" = input_folder_C
)

# ------------------- Output File Paths -------------------
metrics_file   <- file.path(output_folder_files, "results_metrics.csv")
predicted_file <- file.path(output_folder_files, "results_predicted.csv")
cpg_file       <- file.path(output_folder_files, "results_cpg.csv")

# remove old files
for (f in c(metrics_file, predicted_file, cpg_file)) {
  if (file.exists(f)) file.remove(f)
}

# ------------------- Main Loop -------------------
for (input_name in names(input_list)) {
  
  input_folder <- input_list[[input_name]]
  message("STARTING: ", input_name, " -> ", input_folder)
  
  rds_files <- list.files(path = input_folder, pattern = ".rds$", recursive = TRUE, full.names = TRUE)
  n_files <- length(rds_files)
  model_counter <- 1
  
  pb <- txtProgressBar(min = 0, max = n_files, style = 3)
  
  for (i in seq_along(rds_files)) {
    file <- rds_files[i]
    message("Processing: ", file)
    
    tryCatch({
      
      obj <- readRDS(file)
      
      # ------------------- Parse directory structure -------------------
      rel_path <- sub(paste0("^", input_folder, "/"), "", file)
      path_parts <- strsplit(rel_path, "/")[[1]]
      file_stem <- tools::file_path_sans_ext(basename(file))
      
      if (input_name == "P") {
        model <- path_parts[1]
        organ <- path_parts[2]
        rsq_threshold <- path_parts[3]
        filter <- sub(".*_matrix_(.*)_clock_modelled$", "\\1", file_stem)
        tile_methode <- "Point"
        epse <- 1
        cpg_grouping <- "Point"
        
      } else if (input_name == "T") {
        model <- path_parts[1]
        tile_methode <- path_parts[2]
        organ <- path_parts[3]
        rsq_threshold <- path_parts[4]
        epse <- as.numeric(sub(".*Clock_(\\d+)bp_.*", "\\1", file_stem))
        filter <- sub(".*_matrix_(.*)_clock_modelled$", "\\1", file_stem)
        cpg_grouping <- "Tile"
        
      } else if (input_name == "C") {
        model <- path_parts[1]
        organ <- path_parts[2]
        epse <- path_parts[3]
        rsq_threshold <- path_parts[4]
        tile_methode <- "Cluster"
        
        filter <- ifelse(
          grepl("min[0-9]+", file_stem),
          sub(".*_(min[0-9]+)(\\.|_).*", "\\1", file_stem),
          "nofilter"
        )
        cpg_grouping <- "Cluster"
      }
      
      # ------------------- Extract Model Data -------------------
      if (model == "Lasso") {
        
        Train_list <- list(
          data = obj$best_metrics$best_model$metadata_train,
          predicted_ages = obj$best_metrics$best_model$train_predicted_ages,
          mae = obj$best_metrics$best_model$train_mae,
          r_squared = obj$best_metrics$best_model$train_r_squared,
          nb_cpg = obj$final_clock$n_cpgs
        )
        
        Validation_list <- list(
          data = obj$final_validation_metadata,
          predicted_ages = obj$best_metrics$final_val_predicted_ages,
          mae = obj$best_metrics$best_val_mae,
          r_squared = obj$best_metrics$best_val_r_squared
        )
        
      } else {
        
        Train_list <- list(
          predicted_ages = obj$best_metrics$best_model$train_predicted_ages,
          mae = obj$best_metrics$best_model$train_mae,
          r_squared = obj$best_metrics$best_model$train_r2,
          nb_cpg = obj$final_clock$n_cpgs
        )
        
        Validation_list <- list(
          predicted_ages = obj$best_metrics$final_val_predicted_ages,
          mae = obj$best_metrics$best_val_mae,
          r_squared = obj$best_metrics$best_val_r_squared
        )
        
        Train_list$data <- metadata %>% filter(Sanger_new_ID %in% rownames(Train_list$predicted_ages))
        Validation_list$data <- metadata %>% filter(Sanger_new_ID %in% rownames(Validation_list$predicted_ages))
      }
      
      # ------------------- Flatten Function -------------------
      flatten_list <- function(lst, dataset_label, model_id) {
        
        metrics_df <- data.frame(
          cpg_grouping = cpg_grouping,
          tile_methode = tile_methode,
          model_id = model_id,
          dataset = dataset_label,
          model = model,
          organ = organ,
          rsq_threshold = rsq_threshold,
          filter = filter,
          epse = epse,
          mae = lst$mae,
          r_squared_value = as.numeric(lst$r_squared),
          stringsAsFactors = FALSE
        )
        
        if (model == "Lasso") {
          predicted_df <- cbind(
            data.frame(
              cpg_grouping = cpg_grouping,
              tile_methode = tile_methode,
              model_id = model_id,
              dataset = dataset_label,
              model = model,
              organ = organ,
              rsq_threshold = rsq_threshold,
              filter = filter,
              epse = epse,
              stringsAsFactors = FALSE
            ),
            predicted_age = as.numeric(lst$predicted_ages)
          )
        } else {
          predicted_df <- cbind(
            data.frame(
              cpg_grouping = cpg_grouping,
              tile_methode = tile_methode,
              model_id = model_id,
              dataset = dataset_label,
              model = model,
              organ = organ,
              rsq_threshold = rsq_threshold,
              filter = filter,
              epse = epse,
              stringsAsFactors = FALSE
            ),
            predicted_age = as.numeric(lst$predicted_ages[, 1])
          )
        }
        
        cpg_df <- data.frame(
          cpg_grouping = cpg_grouping,
          tile_methode = tile_methode,
          model_id = model_id,
          dataset = dataset_label,
          model = model,
          organ = organ,
          rsq_threshold = rsq_threshold,
          filter = filter,
          epse = epse,
          num_selected_cpg = Train_list$nb_cpg,
          num_candidate_cpg = obj$best_metrics$best_model$num_candidate_cpgs,
          selected_cpgs = paste(obj$best_metrics$best_model$selected_Cpg, collapse = ";"),
          stringsAsFactors = FALSE
        )
        
        return(list(metrics = metrics_df, predicted = predicted_df, cpg = cpg_df))
      }
      
      # ------------------- Flatten Outputs -------------------
      train_flat <- flatten_list(Train_list, "Train", model_counter)
      validation_flat <- flatten_list(Validation_list, "Validation", model_counter)
      
      # ------------------- Plot Paths -------------------
      plot_file_name <- paste0(tools::file_path_sans_ext(basename(rel_path)), "_plot.pdf")
      sub_dirs <- dirname(rel_path)
      plot_output_dir <- file.path(output_folder_plots, cpg_grouping, sub_dirs)
      dir.create(plot_output_dir, recursive = TRUE, showWarnings = FALSE)
      plot_full_path <- file.path(plot_output_dir, plot_file_name)
      
      # ------------------- ggplot2 Combined Plot -------------------
      p_train <- plot_predicted_vs_chron_gg(Train_list, "Train")
      p_val   <- plot_predicted_vs_chron_gg(Validation_list, "Validation")
      
      p_combined <- p_train + p_val + plot_layout(ncol = 2)
      
      ggsave(
        filename = plot_full_path,
        plot = p_combined,
        width = 10,
        height = 5,
        device = "pdf"
      )
      
      # ------------------- Save Metrics -------------------
      append_or_write(train_flat$metrics, metrics_file)
      append_or_write(validation_flat$metrics, metrics_file)
      
      append_or_write(train_flat$predicted, predicted_file)
      append_or_write(validation_flat$predicted, predicted_file)
      
      append_or_write(train_flat$cpg, cpg_file)
      append_or_write(validation_flat$cpg, cpg_file)
      
      model_counter <- model_counter + 1
      gc()
      
    }, error = function(e) {
      message("Skipping file: ", file, " - ", e$message)
    })
    
    setTxtProgressBar(pb, i)
  }
  
  close(pb)
  message("Completed: ", input_name)
}

# ------------------- Load Final Collated Results -------------------
metrics_results   <- read.csv(metrics_file)
predicted_results <- read.csv(predicted_file)
cpg_results       <- read.csv(cpg_file)
