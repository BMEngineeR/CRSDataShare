#!/usr/bin/env Rscript

# Script to create example visualizations for CRS data
# This script generates plots to demonstrate the structure and content of the datasets

# Load required libraries
library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(readr)
library(RColorBrewer)

# Function to create UMAP visualization from Seurat object
plot_umap <- function(seurat_obj, group_by, title, filename) {
  p <- DimPlot(seurat_obj, reduction = "umap", group.by = group_by, label = TRUE, repel = TRUE) +
    ggtitle(title) +
    theme(legend.position = "right") +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16))
  
  # Save the plot
  ggsave(filename, plot = p, width = 12, height = 10, dpi = 300)
  return(p)
}

# Function to create gene expression violin plot
plot_gene_expression <- function(seurat_obj, gene_list, group_by, title, filename) {
  # Ensure genes exist in the dataset
  existing_genes <- intersect(gene_list, rownames(seurat_obj))
  
  if (length(existing_genes) > 0) {
    p <- VlnPlot(seurat_obj, features = existing_genes, group.by = group_by, 
                 pt.size = 0, ncol = min(3, length(existing_genes))) +
      ggtitle(title) +
      theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16))
    
    # Save the plot
    ggsave(filename, plot = p, width = min(15, 5*length(existing_genes)), height = 8, dpi = 300)
    return(p)
  } else {
    cat("None of the specified genes were found in the dataset.\n")
    return(NULL)
  }
}

# Function to visualize cell type proportions by condition
plot_cell_proportions <- function(seurat_obj, cell_type_col, condition_col, title, filename) {
  # Extract metadata
  meta <- seurat_obj@meta.data
  
  # Create contingency table
  cell_counts <- table(meta[[cell_type_col]], meta[[condition_col]])
  
  # Convert to proportions within each condition
  cell_props <- prop.table(cell_counts, margin = 2) * 100
  
  # Convert to data frame for ggplot
  cell_props_df <- as.data.frame(cell_props)
  colnames(cell_props_df) <- c("CellType", "Condition", "Percentage")
  
  # Create the plot
  p <- ggplot(cell_props_df, aes(x = Condition, y = Percentage, fill = CellType)) +
    geom_bar(stat = "identity", position = "fill", width = 0.7) +
    scale_y_continuous(labels = scales::percent_format()) +
    labs(title = title, y = "Percentage of Cells", x = "Condition") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
          legend.position = "right")
  
  # Save the plot
  ggsave(filename, plot = p, width = 10, height = 8, dpi = 300)
  return(p)
}

# Function to create PCA plot showing batch effect
plot_pca_batch <- function(seurat_obj, batch_col, title, filename) {
  p <- DimPlot(seurat_obj, reduction = "pca", group.by = batch_col, 
               label = FALSE, repel = TRUE) +
    ggtitle(title) +
    theme(legend.position = "right") +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16))
  
  # Save the plot
  ggsave(filename, plot = p, width = 10, height = 8, dpi = 300)
  return(p)
}

# Function to analyze GeoMX data and create heatmap of top variable genes
plot_geomx_heatmap <- function(expr_file, anno_file, title, filename, n_genes = 50) {
  # Read the data
  expr <- read.csv(expr_file, row.names = 1)
  anno <- read.csv(anno_file, row.names = 1)
  
  # Identify top variable genes
  gene_vars <- apply(expr, 1, var)
  top_genes <- names(sort(gene_vars, decreasing = TRUE))[1:n_genes]
  
  # Extract expression matrix for top genes
  expr_top <- expr[top_genes, ]
  
  # Prepare data for heatmap
  expr_mat <- as.matrix(expr_top)
  
  # Simplify annotation for visualization
  if ("TissueType" %in% colnames(anno)) {
    anno_col <- anno[colnames(expr_mat), "TissueType", drop = FALSE]
  } else {
    # Use first column as default
    anno_col <- anno[colnames(expr_mat), 1, drop = FALSE]
  }
  
  # Plot heatmap
  if (require("pheatmap")) {
    pheatmap::pheatmap(expr_mat, scale = "row", 
                      clustering_distance_rows = "correlation",
                      clustering_distance_cols = "correlation",
                      annotation_col = anno_col,
                      main = title,
                      filename = filename,
                      width = 12, height = 10)
  } else {
    cat("The pheatmap package is not available. Installing...\n")
    install.packages("pheatmap")
    if (require("pheatmap")) {
      pheatmap::pheatmap(expr_mat, scale = "row", 
                        clustering_distance_rows = "correlation",
                        clustering_distance_cols = "correlation",
                        annotation_col = anno_col,
                        main = title,
                        filename = filename,
                        width = 12, height = 10)
    } else {
      cat("Failed to install pheatmap package. Heatmap cannot be generated.\n")
    }
  }
}

# Main execution
cat("Creating visualizations for CRS data...\n")

# Create output directory for plots
dir.create("CRS_Data_Visualizations", showWarnings = FALSE)

# 1. Plot Public scRNA-seq data
tryCatch({
  cat("Reading public scRNA-seq data...\n")
  public_seurat <- readRDS("CRS_Data/NatureImmunitydataPublic2022/crs.all_seurat_jan24_24.rds")
  
  # UMAP by disease condition
  cat("Plotting UMAP by disease condition...\n")
  plot_umap(public_seurat, "disease", "Public scRNA-seq Data - Clusters by Disease", 
            "CRS_Data_Visualizations/public_umap_disease.png")
  
  # UMAP by cell type
  cat("Plotting UMAP by cell type...\n")
  plot_umap(public_seurat, "celltype", "Public scRNA-seq Data - Cell Types", 
            "CRS_Data_Visualizations/public_umap_celltype.png")
  
  # PCA plot showing batch effect
  cat("Plotting PCA to visualize batch effect...\n")
  plot_pca_batch(public_seurat, "batch", "Public scRNA-seq Data - Batch Effect", 
               "CRS_Data_Visualizations/public_pca_batch.png")
  
  # Cell proportions by disease
  cat("Plotting cell type proportions by disease...\n")
  plot_cell_proportions(public_seurat, "celltype", "disease", 
                      "Cell Type Distribution by Disease Condition", 
                      "CRS_Data_Visualizations/public_cell_proportions.png")
  
  # Gene expression of marker genes (examples)
  marker_genes <- c("CD3E", "CD4", "CD8A", "IL17A", "FOXP3", "CD19", "CD14", "TPSAB1")
  cat("Plotting expression of marker genes...\n")
  plot_gene_expression(public_seurat, marker_genes, "disease", 
                     "Marker Gene Expression by Disease", 
                     "CRS_Data_Visualizations/public_marker_genes.png")
  
  # Clear the large object
  rm(public_seurat)
  gc()
}, error = function(e) {
  cat("Error processing public scRNA-seq data:", e$message, "\n")
})

# 2. Plot In-house Epithelial scRNA-seq data
tryCatch({
  cat("Reading in-house epithelial scRNA-seq data...\n")
  epi_seurat <- readRDS("CRS_Data/InHouseScRNA/Epi_reannotated.rds")
  
  # UMAP by disease condition
  cat("Plotting UMAP by sample type...\n")
  plot_umap(epi_seurat, "subtype", "Epithelial scRNA-seq Data - Clusters by Sample Type", 
            "CRS_Data_Visualizations/epithelial_umap_subtype.png")
  
  # UMAP by cell type
  cat("Plotting UMAP by cell subtype...\n")
  plot_umap(epi_seurat, "Sub_cluster", "Epithelial scRNA-seq Data - Cell Subtypes", 
            "CRS_Data_Visualizations/epithelial_umap_subcluster.png")
  
  # Cell proportions by disease
  cat("Plotting epithelial cell subtypes by disease condition...\n")
  plot_cell_proportions(epi_seurat, "Sub_cluster", "subtype", 
                      "Epithelial Cell Subtype Distribution by Disease Condition", 
                      "CRS_Data_Visualizations/epithelial_subtype_proportions.png")
  
  # Gene expression of epithelial marker genes (examples)
  epi_marker_genes <- c("KRT5", "KRT14", "FOXJ1", "SCGB1A1", "MUC5AC", "ITPKB")
  cat("Plotting expression of epithelial marker genes...\n")
  plot_gene_expression(epi_seurat, epi_marker_genes, "Sub_cluster", 
                     "Epithelial Marker Gene Expression by Cell Subtype", 
                     "CRS_Data_Visualizations/epithelial_marker_genes.png")
  
  # Clear the large object
  rm(epi_seurat)
  gc()
}, error = function(e) {
  cat("Error processing epithelial scRNA-seq data:", e$message, "\n")
})

# 3. Plot In-house Immune Cell scRNA-seq data
tryCatch({
  cat("Reading in-house immune cell scRNA-seq data...\n")
  immune_seurat <- readRDS("CRS_Data/InHouseScRNA/ImmuneCell_processed.rds")
  
  # UMAP by disease condition
  cat("Plotting UMAP by sample type...\n")
  plot_umap(immune_seurat, "subtype", "Immune Cell scRNA-seq Data - Clusters by Sample Type", 
            "CRS_Data_Visualizations/immune_umap_subtype.png")
  
  # UMAP by cell type
  cat("Plotting UMAP by cell type...\n")
  plot_umap(immune_seurat, "Clusters_name", "Immune Cell scRNA-seq Data - Cell Types", 
            "CRS_Data_Visualizations/immune_umap_celltype.png")
  
  # Cell proportions by disease
  cat("Plotting immune cell type proportions by disease condition...\n")
  plot_cell_proportions(immune_seurat, "Clusters_name", "subtype", 
                      "Immune Cell Type Distribution by Disease Condition", 
                      "CRS_Data_Visualizations/immune_celltype_proportions.png")
  
  # Gene expression of immune marker genes (examples)
  immune_marker_genes <- c("CD3E", "CD4", "CD8A", "FOXP3", "CD19", "CD14", "FCGR3A", "MS4A1")
  cat("Plotting expression of immune marker genes...\n")
  plot_gene_expression(immune_seurat, immune_marker_genes, "Clusters_name", 
                     "Immune Marker Gene Expression by Cell Type", 
                     "CRS_Data_Visualizations/immune_marker_genes.png")
  
  # Clear the large object
  rm(immune_seurat)
  gc()
}, error = function(e) {
  cat("Error processing immune cell scRNA-seq data:", e$message, "\n")
})

# 4. Plot GeoMX data
tryCatch({
  cat("Reading GeoMX data...\n")
  
  # Create heatmap of top variable genes
  cat("Creating heatmap of top variable genes in GeoMX data...\n")
  plot_geomx_heatmap("CRS_Data/GeoMX_sequencing_data/countFile_normalized_and_batch_effect_corrected.csv",
                   "CRS_Data/GeoMX_sequencing_data/sampleAnnoFile.csv",
                   "Top Variable Genes in GeoMX Spatial Data",
                   "CRS_Data_Visualizations/geomx_top_variable_genes.png",
                   n_genes = 50)
  
}, error = function(e) {
  cat("Error processing GeoMX data:", e$message, "\n")
})

cat("\nVisualization creation complete. Results saved in 'CRS_Data_Visualizations' directory.\n")