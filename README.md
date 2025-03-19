# CRS Project Data Description

## Overview

This document describes the data available for the Chronic rhinosinusitis (CRS) project, including single-cell RNA sequencing (scRNA-seq) and Geographical Molecular Profiler (GeoMX) datasets. The purpose of this document is to provide a comprehensive understanding of the data formats, structures, and content to facilitate further analysis.

## Data Types

The CRS project contains three main types of data:

1. **Public scRNA-seq Data**: A publicly available dataset from a Nature Immunity 2022 publication.
2. **In-house scRNA-seq Data**: Two datasets representing epithelial cells and immune cells.
3. **In-house GeoMX Data**: Consisting of both sequencing data and imaging data.

## 1. Public scRNA-seq Dataset

### Data Location
- File: `CRS_Data/NatureImmunitydataPublic2022/crs.all_seurat_jan24_24.rds`
- Format: R Seurat object (.rds)
- Access Method: `readRDS()` function in R

### Description
This is a publicly available scRNA-seq dataset from a Nature Immunity 2022 publication related to Chronic rhinosinusitis (CRS). The dataset has been processed and stored as a Seurat object.

### Detailed Analysis
- **Object class**: Seurat
- **Number of cells**: 91,212
- **Number of features (genes)**: 32,643
- **Available assays**: RNA, SCT
- **Dimensional reductions**: PCA, Harmony, UMAP (all available)
- **Disease groups**: 
  - CRSsNP: 21,605 cells
  - Control: 27,550 cells
  - eCRSwNP: 14,240 cells
  - neCRSwNP: 27,817 cells
- **Cell types**: 24 distinct cell types identified
- **Batch information**: 12 different batches (n00, n05, n06, n07, n10, n11, n12, n94, n95, n96, n97, n99)

The public dataset is well-processed, with normalized data, batch correction using Harmony, and dimensionality reduction already performed. Cell clustering and annotation have also been completed, making this dataset ready for downstream analysis.

## 2. In-house scRNA-seq Datasets

### Dataset 1: Epithelial Cells
- File: `CRS_Data/InHouseScRNA/Epi_reannotated.rds`
- Format: R Seurat object (.rds)
- Access Method: `readRDS()` function in R

#### Analysis
- **Object class**: Seurat
- **Number of cells**: 21,833
- **Number of features (genes)**: 24,407
- **Available assays**: RNA, ADT, integrated
- **Dimensional reductions**: PCA, UMAP (both available)
- **Sample types**:
  - Control_Eth: 7,034 cells
  - CRSsNP_Eth: 5,841 cells
  - CRSwNP_Eth: 1,457 cells
  - CRSwNP_NP: 7,501 cells
- **Cell subtypes**: 13 epithelial cell subtypes, including:
  - Basal cells: 2,242 cells
  - Ciliated cells: 4,709 cells
  - Cycling basal cells: 342 cells
  - Deuterosomal cells: 94 cells
  - FOXJ1low Ciliated cells: 812 cells
  - Goblet cells: 1,285 cells
  - Intermediate ciliated cells: 552 cells
  - Ionocytes: 383 cells
  - Mucous cell: 156 cells
  - Secretory cells: 6,544 cells
  - Serous cells: 2,187 cells
  - Suprabasal cells: 2,412 cells
  - Tuft cells: 115 cells

### Dataset 2: Immune Cells
- File: `CRS_Data/InHouseScRNA/ImmuneCell_processed.rds`
- Format: R Seurat object (.rds)
- Access Method: `readRDS()` function in R

#### Analysis
- **Object class**: Seurat
- **Number of cells**: 32,775
- **Number of features (genes)**: 24,407
- **Available assays**: RNA, ADT, integrated
- **Dimensional reductions**: PCA, UMAP (both available)
- **Sample types**:
  - Control_Eth: 6,633 cells
  - CRSsNP_Eth: 10,193 cells
  - CRSwNP_Eth: 4,510 cells
  - CRSwNP_NP: 11,439 cells
- **Cell types**: 13 immune cell types, including:
  - B cells: 2,545 cells
  - CD8 T cells: 6,017 cells
  - Macrophages: 1,502 cells
  - Mast cells: 985 cells
  - Mono/Macrophages: 2,890 cells
  - Naive T cells: 3,002 cells
  - Neutrophils: 2,154 cells
  - NK cells: 1,351 cells
  - NKT cells: 1,049 cells
  - pDCs: 832 cells
  - Plasma cells: 1,439 cells
  - Regulatory T cells: 1,467 cells
  - T cells: 5,540 cells

Both in-house scRNA-seq datasets are well-processed with integrated analysis, dimensionality reduction, and cell type annotation. The datasets include both RNA and Antibody-Derived Tag (ADT) assays, suggesting multi-modal data acquisition that combines transcriptomic and protein-level information.

## 3. GeoMX Datasets

### 3.1 GeoMX Sequencing Data

#### Data Files
- Raw Counts: `CRS_Data/GeoMX_sequencing_data/countFile_raw.csv`
- Processed Counts: `CRS_Data/GeoMX_sequencing_data/countFile_normalized_and_batch_effect_corrected.csv`
- Feature Annotations: `CRS_Data/GeoMX_sequencing_data/featureAnnoFile.csv`
- Sample Annotations: `CRS_Data/GeoMX_sequencing_data/sampleAnnoFile.csv`

#### Analysis
- **Raw counts matrix**: 18,677 genes × 545 samples
- **Processed counts matrix**: 18,676 genes × 537 samples
- **Feature annotation**: 18,815 genes × 10 attributes (including ProbeName, ProbeDisplayName, TargetName, HUGOSymbol, Accessions, GenomeBuild, GenomicPosition, etc.)
- **Sample annotation**: 545 samples × 37 attributes

#### Sample Metadata Structure
The sample metadata includes detailed information about each Area of Interest (AOI):
- **Slide and ROI information**: SlideName, ScanLabel, ROILabel, SegmentLabel
- **Tissue classification**: TissueType (CRSwNP, Spleen, LymphNode, Tonsil, CRSsNP), NewType (neCRSwNP, eCRSwNP, etc.)
- **Anatomical origin**: PathoType (NP, Spleen, LymphNode, Tonsil, UNC)
- **AOI metrics**: AOISurfaceArea, AOINucleiCount
- **ROI coordinates**: ROICoordinateX, ROICoordinateY
- **Sequencing metrics**: RawReads, AlignedReads, DeduplicatedReads, TrimmedReads, StitchedReads, SequencingSaturation

Each sample appears to represent a specific segment (MAC, IMM, EPI) within a tissue microarray (TMA) spot, suggesting that the data captures spatial information about different cellular compartments within the tissue.

### 3.2 GeoMX Imaging Data

#### Data Location
- Directory: `CRS_Data/GeoMX_image/`
- Number of TMA image directories: 183
- Sample TMA directories: TMA016001, TMA016002, TMA016003, etc.

#### Structure of Image Data
Each TMA directory contains a set of TIFF images representing different staining channels and segmentation information:

**Staining Channels**:
- CD45.tiff: Staining for CD45, a common leukocyte marker
- CD68.tiff: Staining for CD68, a macrophage marker
- PanCK.tiff: Pan-cytokeratin staining for epithelial cells
- SYTO13.tiff: Nuclear staining
- membrane.tiff: Membrane staining
- nuclear.tiff: Alternative nuclear staining

**Segmentation Files**:
- MESMER_mask.tiff: Cell segmentation mask generated by MESMER algorithm
- seg_outline.tiff: Segmentation outlines
- seg_overlay.tiff: Segmentation overlay on the original image

This imaging data provides the spatial context for the sequencing data, allowing for correlation between gene expression patterns and tissue morphology. The segmentation files indicate that cell-level resolution is available for spatial analysis.

## Integration Potential

The CRS project datasets offer significant opportunities for integrated analysis:

1. **Public vs. In-house scRNA-seq Comparison**: The public dataset (91,212 cells) can be compared with the in-house datasets (21,833 epithelial cells + 32,775 immune cells) to validate findings and identify novel cell populations or disease-specific signatures.

2. **Cell Type-Specific Analysis**: The detailed cell type annotations in both scRNA-seq datasets allow for focused analysis of specific cell populations across conditions (Control, CRSsNP, CRSwNP).

3. **Spatial-Molecular Integration**: The GeoMX data provides spatial context that complements the single-cell resolution of the scRNA-seq data. Gene expression patterns identified in scRNA-seq can be mapped onto the spatial data to understand their tissue context.

4. **Tissue Microenvironment Analysis**: The segmentation of GeoMX data into different compartments (MAC, IMM, EPI) enables analysis of cell-cell interactions and tissue microenvironment effects on gene expression.

## Recommended Analysis Workflow

Based on the data structure and content, the following analysis workflow is recommended:

1. **Quality Assessment**:
   - Verify cell type annotations across datasets
   - Check for batch effects and integration quality
   - Assess coverage of genes of interest

2. **Differential Expression Analysis**:
   - Compare disease conditions within each dataset
   - Identify cell type-specific markers
   - Explore pathway enrichment

3. **Spatial Analysis**:
   - Map cell types from scRNA-seq to spatial regions in GeoMX
   - Analyze compartment-specific gene expression
   - Investigate spatial patterns of disease-associated genes

4. **Integrative Analysis**:
   - Correlate findings between public and in-house datasets
   - Develop predictive models using multiple data modalities
   - Identify key regulators and potential therapeutic targets

## Potential Research Questions

The CRS datasets are particularly well-suited to address several key research questions:

1. **Cellular Heterogeneity in CRS**:
   - How do cell type compositions differ between healthy controls and CRS patients?
   - Are there novel cell states or subtypes unique to CRS conditions?
   - How does the epithelial cell landscape differ between CRSsNP and CRSwNP?

2. **Inflammatory Mechanisms**:
   - What are the key inflammatory pathways activated in different CRS subtypes?
   - Which immune cell populations drive eosinophilic vs. non-eosinophilic CRSwNP?
   - What are the signaling interactions between epithelial and immune cells in the nasal mucosa?

3. **Spatial Organization**:
   - How are different cell types spatially organized in CRS tissues?
   - Do specific immune cell populations aggregate in particular tissue regions?
   - How does the spatial architecture of the tissue relate to disease severity?

4. **Biomarker Discovery**:
   - Can we identify gene signatures that distinguish CRS subtypes?
   - Are there potential therapeutic targets expressed in specific cell populations?
   - Can spatial gene expression patterns predict disease progression or treatment response?

5. **Comparison with Literature**:
   - How do our findings compare with published CRS studies?
   - Can we validate previously reported disease mechanisms?
   - What novel insights does our multi-modal approach provide?

## Example Analysis Code

Below are code snippets demonstrating key analyses that can be performed with these datasets:

### 1. Loading and Exploring scRNA-seq Data

```R
# Load Seurat object
public_seurat <- readRDS("CRS_Data/NatureImmunitydataPublic2022/crs.all_seurat_jan24_24.rds")

# Basic exploration
print(public_seurat)
table(public_seurat$disease)
table(public_seurat$celltype)

# Visualize UMAP
DimPlot(public_seurat, reduction = "umap", group.by = "disease", label = TRUE)
DimPlot(public_seurat, reduction = "umap", group.by = "celltype", label = TRUE)

# Cell type proportions by disease
prop.table(table(public_seurat$celltype, public_seurat$disease), margin = 2) * 100
```

### 2. Differential Expression Analysis

```R
# Differential expression between conditions for a specific cell type
cell_type <- "Secretory cells"
Idents(public_seurat) <- public_seurat$celltype
markers <- FindMarkers(public_seurat, ident.1 = "CRSwNP", ident.2 = "Control", 
                        subset.ident = cell_type, 
                        min.pct = 0.25, logfc.threshold = 0.25)
head(markers, n = 20)
```

### 3. GeoMX Data Analysis

```R
# Load GeoMX data
raw_counts <- read.csv("CRS_Data/GeoMX_sequencing_data/countFile_raw.csv")
processed_counts <- read.csv("CRS_Data/GeoMX_sequencing_data/countFile_normalized_and_batch_effect_corrected.csv")
feature_anno <- read.csv("CRS_Data/GeoMX_sequencing_data/featureAnnoFile.csv")
sample_anno <- read.csv("CRS_Data/GeoMX_sequencing_data/sampleAnnoFile.csv")

# Heatmap of top variable genes
library(pheatmap)
top_var_genes <- names(sort(apply(processed_counts, 1, var), decreasing = TRUE))[1:100]
pheatmap(processed_counts[top_var_genes, ], 
         annotation_col = sample_anno[, c("TissueType", "NewType")],
         show_rownames = FALSE, show_colnames = FALSE,
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation")
```

## Conclusion

The CRS project data represents a comprehensive resource for studying Chronic rhinosinusitis through multiple modalities (scRNA-seq and spatial transcriptomics) and sources (public and in-house). The datasets are well-processed with detailed annotations, making them ready for advanced integrative analyses to understand the cellular and molecular basis of CRS pathology.
