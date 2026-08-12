# ==========================================================================
# PIPELINE STAGE: 01 / Annotation
# ==========================================================================
# Builds the per-CpG methylation matrices (Point + Tile CpG-island clocks)
# from Bismark .cov.gz coverage files. Run this first; its outputs feed
# stage 02 (candidate clock construction, steps 1-5).
# 
# Shared verbatim between both aims of the study (clock design benchmarking
# and cross-tissue ageing biology) - both use the exact same underlying
# methylation matrices.
# ==========================================================================

# Set the root directory of this project (edit here, or export EPICLOCK_ROOT
# before launching R / sbatch). All absolute paths below are built from this.
PROJECT_ROOT <- Sys.getenv("EPICLOCK_ROOT", unset = normalizePath("."))

# --- R Script: Generate per-CpG methylation matrix for TILE and POINT epigenetic clock construction ---

# This script processes Bismark .cov.gz files (bismarkCoverage pipeline),
# extracts per-CpG methylation values, applies coverage filtering,
# unites samples retaining sites present in ≥90% samples,
# and outputs methylation percentage matrices filtered at multiple stringency levels.

# Load required packages ----

### methylkit requires a specific verison of data.table 
#install.packages("https://cran.r-project.org/src/contrib/Archive/data.table/data.table_1.14.2.tar.gz", repos = NULL, type = "source")

# --- Dual Clock Methylation Matrix Generator (Point + Tile CpG Island Clocks) ---

# had some issues with temporary memory so forcing a directory
# Redirect R's temporary directory to avoid /tmp issues
Sys.setenv(TMPDIR = file.path(PROJECT_ROOT, "../tmp"))


# Ensure BiocManager is installed
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", repos = "https://cran.rstudio.com/")

bioc_packages <- c(
  "methylKit",
  "BSgenome.Drerio.UCSC.danRer11",
  "GenomicRanges",
  "IRanges"
)

for (pkg in bioc_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg, ask = FALSE)
  }
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}


# CRAN packages to check/install
cran_packages <- c(
  "data.table", "dplyr", "readr"
)

invisible(lapply(cran_packages, function(pkg) {
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}))

# --- Configuration ----

output_dir <- "Results/Annotations"
base_dir <-  "data/cov_files"
cpg_island_bed <- "data/danRer11_CpGislands.bed"
# Create output directory if missing
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)


tile_sizes <- c(50, 100, 500, 1000, 2000, 5000, 25000)
filter_thresholds <- list("nofilter" = NULL, "min5" = 5, "min10" = 10, "min20" = 20)

# --- Load CpG islands as GRanges ----
load_cpg_islands <- function(bed_file) {
  df <- readr::read_tsv(bed_file, col_names = c("chr", "start", "end"))
  df$start <- df$start + 1  # BED format is 0-based
  GRanges(seqnames = df$chr, ranges = IRanges(start = df$start, end = df$end))
}
cpg_islands <- load_cpg_islands(cpg_island_bed)

# --- Per-folder processing ----
process_cov_folder <- function(folder_path, output_base) {
  message("Processing folder: ", folder_path)
  
  cov_files <- list.files(folder_path, pattern = "\\.cov\\.gz$", full.names = TRUE)
  if (length(cov_files) == 0) {
    message("No .cov.gz files found in ", folder_path)
    return(NULL)
  }
  
  sample_ids <- sub("\\.sorted.*", "", basename(cov_files))
  treatment_vector <- rep(0, length(cov_files))
  
  meth_raw_list <- lapply(seq_along(cov_files), function(i) {
    methRead(
      location = cov_files[i],
      sample.id = sample_ids[i],
      assembly = "danRer11",
      treatment = treatment_vector[i],
      context = "CpG",
      pipeline = "bismarkCoverage",
      mincov = 5
    )
  })
  meth_raw <- new("methylRawList", meth_raw_list)
  
  n_samples <- length(sample_ids)
  min_samples <- as.integer(ceiling(n_samples * 0.9))
  
  coverage_filters <- list()
  
  for (label in names(filter_thresholds)) {
    thresh <- filter_thresholds[[label]]
    filtered_raw <- if (is.null(thresh)) {
      meth_raw
    } else {
      new("methylRawList", lapply(meth_raw, function(x) {
        filterByCoverage(x, lo.count = thresh, hi.perc = 99.9)
      }))
    }
    united <- methylKit::unite(filtered_raw, destrand = FALSE, min.per.group = min_samples)
    coverage_filters[[label]] <- united
  }
  
  # Prepare output subdirectories
  #Extract organ name
  organ_name <- basename(folder_path)
  
  # Prepare output subdirectories with flipped hierarchy
  point_dir       <- file.path(output_base, "point_clock", organ_name)
  cpg_tile_dir    <- file.path(output_base, "cpg_island_tile_clock", organ_name)
  genome_tile_dir <- file.path(output_base, "whole_genome_tile_clock", organ_name)
  
  dir.create(point_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(cpg_tile_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(genome_tile_dir, recursive = TRUE, showWarnings = FALSE)
  
  # --- Point Clock ---
  lapply(names(coverage_filters), function(label) {
    meth_obj <- coverage_filters[[label]]
    df <- getData(meth_obj)
    
    coords <- df %>% dplyr::select(chr, start, end, strand)
    sample_ids <- meth_obj@sample.ids
    
    perc_meth_mat <- matrix(NA_real_, nrow = nrow(df), ncol = length(sample_ids),
                            dimnames = list(NULL, sample_ids))
    
    for (i in seq_along(sample_ids)) {
      coverage_col <- paste0("coverage", i)
      numCs_col <- paste0("numCs", i)
      perc_meth_mat[, i] <- ifelse(df[[coverage_col]] > 0,
                                   (df[[numCs_col]] / df[[coverage_col]]) * 100,
                                   NA_real_)
    }
    
    final_df <- cbind(coords, as.data.frame(perc_meth_mat))
    out_file <- file.path(point_dir, paste0("pointClock_methylation_matrix_", label, ".tsv"))
    
    write.table(final_df, file = out_file, sep = "\t", quote = FALSE, row.names = FALSE)
    message("Wrote point clock matrix: ", out_file)
    
    gc()
  })
  
  # --- CpG Island Tile Clock ---
  for (tile_width in tile_sizes) {
    message("Tiling CpG islands: ", tile_width, " bp")
    tiled_islands <- unlist(slidingWindows(cpg_islands, width = tile_width, step = tile_width))
    
    for (label in names(coverage_filters)) {
      meth_obj <- coverage_filters[[label]]
      meth_obj$chr <- paste0("chr", meth_obj$chr)  # Ensure chr format matches
      
      regional <- regionCounts(meth_obj, tiled_islands)
      df <- getData(regional)
      
      coords <- df %>% dplyr::select(chr, start, end, strand)
      sample_ids <- regional@sample.ids
      
      perc_meth_mat <- matrix(NA_real_, nrow = nrow(df), ncol = length(sample_ids),
                              dimnames = list(NULL, sample_ids))
      
      for (i in seq_along(sample_ids)) {
        coverage_col <- paste0("coverage", i)
        numCs_col <- paste0("numCs", i)
        perc_meth_mat[, i] <- ifelse(df[[coverage_col]] > 0,
                                     (df[[numCs_col]] / df[[coverage_col]]) * 100,
                                     NA_real_)
      }
      
      final_df <- cbind(coords, as.data.frame(perc_meth_mat))
      out_file <- file.path(cpg_tile_dir,
                            paste0("tileClock_", tile_width, "bp_methylation_matrix_", label, ".tsv"))
      
      write.table(final_df, file = out_file, sep = "\t", quote = FALSE, row.names = FALSE)
      message("Wrote CpG island tile clock matrix: ", out_file)
      
      gc()
    }
  }
  
  # --- Whole Genome Tile Clock ---
  # Define chromosome lengths: you need to provide this or infer it
  
  # Example: get chromosome lengths from your data or genome reference
  # Here, infer from methylRawList assuming consistent chromosomes
  # Convert chromosome names in meth_raw to "chr" prefix as needed
  meth_obj_for_chr <- coverage_filters[[1]]
  meth_obj_for_chr$chr <- paste0("chr", meth_obj_for_chr$chr)
  # need the chromosome lengths 
  chrom_lengths <- seqlengths(BSgenome.Drerio.UCSC.danRer11)
  # removing all non chromosomal data and only keeping 25 chr data 
  chrom_lengths <- chrom_lengths[1:25]
  
  
  # If seqlengths not set, you must provide chrom lengths as named vector:
  # e.g. chrom_lengths <- c(chr1=XXXX, chr2=YYYY, ...)
  
  for (tile_width in tile_sizes) {
    message("Tiling whole genome: ", tile_width, " bp")
    
    # Create whole genome tiles
    genome_tiles <- tileGenome(chrom_lengths,
                               tilewidth = tile_width,
                               cut.last.tile.in.chrom = TRUE)
    
    for (label in names(coverage_filters)) {
      meth_obj <- coverage_filters[[label]]
      meth_obj$chr <- paste0("chr", meth_obj$chr)  # Ensure consistent chr format
      
      regional <- regionCounts(meth_obj, genome_tiles)
      df <- getData(regional)
      
      coords <- df %>% dplyr::select(chr, start, end, strand)
      sample_ids <- regional@sample.ids
      
      perc_meth_mat <- matrix(NA_real_, nrow = nrow(df), ncol = length(sample_ids),
                              dimnames = list(NULL, sample_ids))
      
      for (i in seq_along(sample_ids)) {
        coverage_col <- paste0("coverage", i)
        numCs_col <- paste0("numCs", i)
        perc_meth_mat[, i] <- ifelse(df[[coverage_col]] > 0,
                                     (df[[numCs_col]] / df[[coverage_col]]) * 100,
                                     NA_real_)
      }
      
      final_df <- cbind(coords, as.data.frame(perc_meth_mat))
      out_file <- file.path(genome_tile_dir,
                            paste0("wholeGenomeTileClock_", tile_width, "bp_methylation_matrix_", label, ".tsv"))
      
      write.table(final_df, file = out_file, sep = "\t", quote = FALSE, row.names = FALSE)
      message("Wrote whole genome tile clock matrix: ", out_file)
      
      gc()
    }
  }
}






# --- Traverse and process all folders ---

all_dirs <- list.dirs(path = base_dir, recursive = FALSE, full.names = TRUE)

# code for testing with smalled folder- brain
#all_dirs1 <- all_dirs[2]

all_dirs1 <- all_dirs

for (folder in all_dirs1) {
  process_cov_folder(folder, output_dir)
}
