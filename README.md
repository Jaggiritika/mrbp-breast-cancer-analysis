# MRBP Breast Cancer Analysis Pipeline

## Overview

This R pipeline performs comprehensive analysis of Mitochondrial RNA Binding Proteins (MRBPs) in breast cancer using The Cancer Genome Atlas (TCGA) BRCA dataset. The analysis includes differential expression analysis, survival analysis, pathway enrichment, and visualization techniques.

## Research Context

Mitochondrial RNA Binding Proteins (MRBPs) play crucial roles in mitochondrial RNA metabolism, processing, and regulation. These proteins are involved in RNA splicing, stability, and translation within mitochondria. This pipeline investigates the expression patterns and clinical significance of MRBPs in breast cancer, providing insights into potential therapeutic targets and biomarkers.

## Features

- **Data Acquisition**: Downloads and processes TCGA-BRCA transcriptome data
- **Differential Expression Analysis**: Identifies significantly differentially expressed MRBPs
- **Survival Analysis**: Correlates MRBP expression with patient survival outcomes
- **Pathway Enrichment**: GO and Reactome pathway analysis
- **Visualization**: Heatmaps, PCA plots, and Venn diagrams
- **Clinical Correlation**: Links gene expression to clinical outcomes

## Prerequisites

### Required R Packages

```r
# Core Bioconductor packages
library(TCGAbiolinks)
library(SummarizedExperiment)
library(DESeq2)
library(org.Hs.eg.db)
library(AnnotationDbi)

# Visualization packages
library(pheatmap)
library(VennDiagram)
library(survminer)
library(survival)

# Pathway analysis packages
library(clusterProfiler)
library(ReactomePA)

# Data manipulation
library(dplyr)
library(grid)
```

### Installation

```r
# Install BiocManager if not already installed
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Install required packages
BiocManager::install(c("TCGAbiolinks", "SummarizedExperiment", "DESeq2", 
                      "org.Hs.eg.db", "AnnotationDbi", "clusterProfiler", 
                      "ReactomePA", "survminer", "pheatmap", "VennDiagram"))
```

## Usage

1. **Data Acquisition**: Downloads TCGA-BRCA transcriptome data
2. **Data Processing**: Filters for primary tumor and normal tissue samples
3. **Differential Expression**: Identifies MRBPs with significant expression changes
4. **Survival Analysis**: Analyzes correlation between MRBP expression and survival
5. **Pathway Analysis**: Performs GO and Reactome enrichment analysis
6. **Visualization**: Creates heatmaps, PCA plots, and Venn diagrams

## Output Files

- `DEG_MRBPs.csv`: Differentially expressed MRBPs
- `GO_Enrichment_MRBPs.csv`: GO enrichment results
- `GO_66.csv`: GO analysis for 66 matched MRBPs
- `reactomePA.csv`: Reactome pathway enrichment results
- `MRBP_DE_results.csv`: Complete MRBP differential expression results
- `matched_gene_names.csv`: Gene symbols for matched MRBPs

## Key Findings

The pipeline identifies 66 MRBPs that are differentially expressed in breast cancer compared to normal tissue. These genes show significant enrichment in mitochondrial and metabolic pathways, suggesting their potential role in cancer metabolism and progression.

## Citation

If you use this pipeline in your research, please cite:

- TCGAbiolinks: Colaprico A, et al. (2016) TCGAbiolinks: an R/Bioconductor package for integrative analysis of TCGA data. Nucleic Acids Research, 44(8):e71.
- DESeq2: Love MI, et al. (2014) Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biology, 15:550.



## Contact

For questions or issues, please open an issue on GitHub or contact the maintainer.

## Acknowledgments

- TCGA Research Network for providing the breast cancer data
- Bioconductor community for the excellent bioinformatics packages
- R community for the robust statistical computing environment 