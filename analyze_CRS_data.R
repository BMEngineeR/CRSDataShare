#!/usr/bin/env Rscript

# Script to analyze CRS project data
# Author: [Your Name]
# Date: [Current Date]

# Load required libraries
library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(readr)

# Set working directory - already set to the shared_data directory
# setwd("/bmbl_data/yuzhou/collaborative/Sizun_lab/CRS/shared_data")

# Function to analyze Seurat objects
analyze_seurat_object <- function(seurat_obj, dataset_name) {
  cat("\n\n## Analysis of", dataset_name, "\n\n")
  
  # Basic object information
  cat("- Object class:", class(seurat_obj)[1], "\n")
  cat("- Number of cells:", ncol(seurat_obj), "\n")
  cat("- Number of features:", nrow(seurat_obj), "\n")
  
  # Check available assays
  cat("- Available assays:", paste(names(seurat_obj@assays), collapse = ", "), "\n")
  
  # Check if dimensional reductions are available
  dimred <- names(seurat_obj@reductions)
  cat("- Available dimensional reductions:", ifelse(length(dimred) > 0, paste(dimred, collapse = ", "), "None"), "\n")
  
  # Summarize metadata
  meta_cols <- colnames(seurat_obj@meta.data)
  cat("- Metadata columns:", paste(meta_cols, collapse = ", "), "\n")
  
  # Show a sample of metadata for categorical columns
  cat("- Summary of categorical metadata:\n")
  for (col in meta_cols) {
    if (is.factor(seurat_obj@meta.data[[col]]) || is.character(seurat_obj@meta.data[[col]])) {
      counts <- table(seurat_obj@meta.data[[col]])
      if (length(counts) < 20) {  # Only show if not too many categories
        cat("  - ", col, ":", "\n")
        for (cat_name in names(counts)) {
          cat("    - ", cat_name, ": ", counts[cat_name], " cells\n", sep="")
        }
      } else {
        cat("  - ", col, ": [Too many categories to display, total: ", length(counts), "]\n", sep="")
      }
    }
  }
  
  # Return key metrics in a list for later use
  return(list(
    dataset = dataset_name,
    cell_count = ncol(seurat_obj),
    feature_count = nrow(seurat_obj),
    assays = names(seurat_obj@assays),
    reductions = dimred,
    metadata_cols = meta_cols
  ))
}

# Function to analyze GeoMX data
analyze_geomx_data <- function() {
  cat("\n\n## Analysis of GeoMX Sequencing Data\n\n")
  
  # Read sequencing data without using row.names to avoid duplicate issues
  raw_counts <- read.csv("CRS_Data/GeoMX_sequencing_data/countFile_raw.csv", row.names = NULL)
  processed_counts <- read.csv("CRS_Data/GeoMX_sequencing_data/countFile_normalized_and_batch_effect_corrected.csv", row.names = NULL)
  feature_anno <- read.csv("CRS_Data/GeoMX_sequencing_data/featureAnnoFile.csv", row.names = NULL)
  sample_anno <- read.csv("CRS_Data/GeoMX_sequencing_data/sampleAnnoFile.csv", row.names = NULL)
  
  # Extract gene IDs to the first column (if it's already there)
  gene_ids_col <- colnames(raw_counts)[1]
  
  # Basic stats for counts
  cat("- Raw counts matrix dimensions:", nrow(raw_counts), "genes ×", ncol(raw_counts) - 1, "samples\n")
  cat("- Processed counts matrix dimensions:", nrow(processed_counts), "genes ×", ncol(processed_counts) - 1, "samples\n")
  cat("- Feature annotation dimensions:", nrow(feature_anno), "genes ×", ncol(feature_anno) - 1, "attributes\n")
  cat("- Sample annotation dimensions:", nrow(sample_anno), "samples ×", ncol(sample_anno) - 1, "attributes\n")
  
  # Analyze sample annotations - first column in feature_anno is likely gene IDs
  cat("\n### Feature Annotation Summary\n\n")
  
  # Get first few rows of feature annotations to show structure
  cat("- First 5 rows of feature annotation data:\n")
  feature_head <- head(feature_anno, 5)
  for (i in 1:nrow(feature_head)) {
    cat("  Row", i, ":", paste(names(feature_head), "=", as.character(feature_head[i,]), collapse = ", "), "\n")
  }
  
  # Analyze sample annotations
  cat("\n### Sample Metadata Summary\n\n")
  
  # Summary of metadata columns
  sample_cols <- colnames(sample_anno)
  cat("- Sample metadata columns:", paste(sample_cols, collapse = ", "), "\n\n")
  
  # Get first few rows of sample annotations to show structure
  cat("- First 5 rows of sample annotation data:\n")
  sample_head <- head(sample_anno, 5)
  for (i in 1:nrow(sample_head)) {
    cat("  Row", i, ":", paste(names(sample_head), "=", as.character(sample_head[i,]), collapse = ", "), "\n")
  }
  
  # Look for patient and TMA information in column names
  patient_col <- grep("patient", sample_cols, ignore.case = TRUE, value = TRUE)
  tma_col <- grep("tma", sample_cols, ignore.case = TRUE, value = TRUE)
  
  # Count unique patients and TMAs if columns exist
  if (length(patient_col) > 0) {
    cat("\n- Number of unique patients:", length(unique(sample_anno[[patient_col[1]]])), "\n")
    cat("- Patient IDs:", paste(head(unique(sample_anno[[patient_col[1]]]), 10), collapse = ", "), 
        ifelse(length(unique(sample_anno[[patient_col[1]]])) > 10, "...", ""), "\n")
  }
  
  if (length(tma_col) > 0) {
    cat("\n- Number of unique TMAs:", length(unique(sample_anno[[tma_col[1]]])), "\n")
    cat("- TMA IDs:", paste(head(unique(sample_anno[[tma_col[1]]]), 10), collapse = ", "), 
        ifelse(length(unique(sample_anno[[tma_col[1]]])) > 10, "...", ""), "\n")
  }
  
  # If specific columns are not found, try to infer from data
  if (length(patient_col) == 0 || length(tma_col) == 0) {
    cat("\n- Note: Explicit patient or TMA columns not found. Searching for relevant columns...\n")
    
    # Check all character/factor columns for potential patient/TMA info
    for (col in sample_cols) {
      if (is.character(sample_anno[[col]]) || is.factor(sample_anno[[col]])) {
        unique_vals <- unique(sample_anno[[col]])
        if (length(unique_vals) > 1 && length(unique_vals) < 100) {  # Reasonable number for patients/TMAs
          # Check if values look like patient IDs or TMA IDs
          if (any(grepl("TMA", unique_vals, ignore.case = TRUE)) || 
              any(grepl("tma", unique_vals, ignore.case = TRUE))) {
            cat("- Column", col, "may contain TMA information.\n")
            cat("  Values (sample):", paste(head(unique_vals, 5), collapse = ", "), "...\n")
          } else if (any(grepl("P", unique_vals)) || 
                    any(grepl("patient", unique_vals, ignore.case = TRUE))) {
            cat("- Column", col, "may contain patient information.\n")
            cat("  Values (sample):", paste(head(unique_vals, 5), collapse = ", "), "...\n")
          }
        }
      }
    }
  }
  
  return(list(
    raw_count_dims = c(nrow(raw_counts), ncol(raw_counts)),
    processed_count_dims = c(nrow(processed_counts), ncol(processed_counts)),
    feature_anno_dims = c(nrow(feature_anno), ncol(feature_anno)),
    sample_anno_dims = c(nrow(sample_anno), ncol(sample_anno)),
    sample_cols = sample_cols
  ))
}

# Function to analyze GeoMX imaging data
analyze_geomx_images <- function() {
  cat("\n\n## Analysis of GeoMX Imaging Data\n\n")
  
  # Get list of TMA directories
  image_dir <- "CRS_Data/GeoMX_image"
  tma_dirs <- list.dirs(image_dir, full.names = FALSE, recursive = FALSE)
  cat("- Number of TMA image directories:", length(tma_dirs), "\n")
  cat("- TMA directories (sample):", paste(head(tma_dirs, 10), collapse = ", "), "...\n\n")
  
  # Choose first TMA as example if available
  if (length(tma_dirs) > 0) {
    example_tma <- tma_dirs[1]
    cat("### Example TMA Structure:", example_tma, "\n\n")
    
    # List files in the example TMA directory
    example_files <- list.files(file.path(image_dir, example_tma), recursive = TRUE, full.names = FALSE)
    cat("- Number of files:", length(example_files), "\n")
    
    # Show file extensions
    file_exts <- table(tools::file_ext(example_files))
    cat("- File types:", "\n")
    for (ext in names(file_exts)) {
      cat("  - ", ext, ": ", file_exts[ext], " files\n", sep="")
    }
    
    # Show directory structure
    example_dirs <- list.dirs(file.path(image_dir, example_tma), full.names = FALSE, recursive = TRUE)
    if (length(example_dirs) > 1) { # Skip the root directory which is ""
      cat("- Directory structure:", "\n")
      for (dir in example_dirs[-1]) { # Skip the first empty directory name
        cat("  - ", dir, "\n", sep="")
      }
    }
    
    # Show file list
    cat("- Files in", example_tma, ":", "\n")
    for (file in example_files) {
      cat("  - ", file, "\n", sep="")
    }
    
    # Describe image content
    cat("\n### Image Content Description\n\n")
    cat("Each TMA directory contains several image files in TIFF format that represent different staining channels and segmentation masks:\n\n")
    cat("- **Staining Channels**:\n")
    cat("  - CD45.tiff: Staining for CD45, a common leukocyte marker\n")
    cat("  - CD68.tiff: Staining for CD68, a macrophage marker\n")
    cat("  - PanCK.tiff: Pan-cytokeratin staining for epithelial cells\n")
    cat("  - SYTO13.tiff: Nuclear staining\n")
    cat("  - membrane.tiff: Membrane staining\n")
    cat("  - nuclear.tiff: Alternative nuclear staining\n\n")
    cat("- **Segmentation Files**:\n")
    cat("  - MESMER_mask.tiff: Cell segmentation mask generated by MESMER algorithm\n")
    cat("  - seg_outline.tiff: Segmentation outlines\n")
    cat("  - seg_overlay.tiff: Segmentation overlay on the original image\n")
  }
  
  return(list(
    tma_count = length(tma_dirs),
    tma_names = tma_dirs
  ))
}

# Main execution
cat("# CRS Data Analysis Results\n\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

# 1. Public scRNA-seq data
cat("# 1. Public scRNA-seq Dataset\n")
tryCatch({
  public_seurat <- readRDS("CRS_Data/NatureImmunitydataPublic2022/crs.all_seurat_jan24_24.rds")
  public_results <- analyze_seurat_object(public_seurat, "Public scRNA-seq (NatureImmunity 2022)")
  # Clear the large object
  rm(public_seurat)
  gc()
}, error = function(e) {
  cat("Error loading public scRNA-seq data:", e$message, "\n")
})

# 2. In-house scRNA-seq data
cat("\n# 2. In-house scRNA-seq Datasets\n")

# Epithelial cells
tryCatch({
  epi_seurat <- readRDS("CRS_Data/InHouseScRNA/Epi_reannotated.rds")
  epi_results <- analyze_seurat_object(epi_seurat, "In-house Epithelial scRNA-seq")
  # Clear the large object
  rm(epi_seurat)
  gc()
}, error = function(e) {
  cat("Error loading epithelial scRNA-seq data:", e$message, "\n")
})

# Immune cells
tryCatch({
  immune_seurat <- readRDS("CRS_Data/InHouseScRNA/ImmuneCell_processed.rds")
  immune_results <- analyze_seurat_object(immune_seurat, "In-house Immune Cell scRNA-seq")
  # Clear the large object
  rm(immune_seurat)
  gc()
}, error = function(e) {
  cat("Error loading immune cell scRNA-seq data:", e$message, "\n")
})

# 3. GeoMX data
cat("\n# 3. GeoMX Datasets\n")

# Sequencing data
tryCatch({
  geomx_seq_results <- analyze_geomx_data()
}, error = function(e) {
  cat("Error analyzing GeoMX sequencing data:", e$message, "\n")
})

# Imaging data
tryCatch({
  geomx_img_results <- analyze_geomx_images()
}, error = function(e) {
  cat("Error analyzing GeoMX imaging data:", e$message, "\n")
})

cat("\n\nAnalysis complete.\n")