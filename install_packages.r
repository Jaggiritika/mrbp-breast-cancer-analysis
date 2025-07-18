# =============================================================================
# Package Installation Script for MRBP Breast Cancer Analysis Pipeline
# =============================================================================
#
# This script installs all required R packages for the MRBP analysis pipeline.
# Run this script before running the main pipeline.
#
# Author: [Your Name]
# Date: [Current Date]
# =============================================================================

# =============================================================================
# 1. Install BiocManager (if not already installed)
# =============================================================================

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
  cat("BiocManager installed successfully.\n")
} else {
  cat("BiocManager is already installed.\n")
}

# =============================================================================
# 2. Install Bioconductor Packages
# =============================================================================

cat("Installing Bioconductor packages...\n")

# Core data analysis packages
BiocManager::install(c(
  "TCGAbiolinks",        # TCGA data access
  "SummarizedExperiment", # Data structure for genomic data
  "DESeq2",              # Differential expression analysis
  "org.Hs.eg.db",        # Human gene annotation database
  "AnnotationDbi"         # Annotation database interface
), update = FALSE)

# Pathway analysis packages
BiocManager::install(c(
  "clusterProfiler",      # Gene set enrichment analysis
  "ReactomePA"           # Reactome pathway analysis
), update = FALSE)

# =============================================================================
# 3. Install CRAN Packages
# =============================================================================

cat("Installing CRAN packages...\n")

# Visualization packages
install.packages(c(
  "pheatmap",            # Heatmap visualization
  "VennDiagram",         # Venn diagram creation
  "survminer",           # Survival analysis visualization
  "survival"             # Survival analysis
))

# Data manipulation packages
install.packages(c(
  "dplyr",               # Data manipulation
  "grid"                 # Grid graphics
))

# =============================================================================
# 4. Install Additional Packages (if needed)
# =============================================================================

# Install biomaRt for gene annotation (if not already installed)
if (!requireNamespace("biomaRt", quietly = TRUE)) {
  BiocManager::install("biomaRt")
  cat("biomaRt installed successfully.\n")
} else {
  cat("biomaRt is already installed.\n")
}

# =============================================================================
# 5. Verify Installation
# =============================================================================

cat("Verifying package installation...\n")

# List of required packages
required_packages <- c(
  "TCGAbiolinks",
  "SummarizedExperiment", 
  "DESeq2",
  "org.Hs.eg.db",
  "AnnotationDbi",
  "clusterProfiler",
  "ReactomePA",
  "pheatmap",
  "VennDiagram",
  "survminer",
  "survival",
  "dplyr",
  "grid",
  "biomaRt"
)

# Check if all packages are available
missing_packages <- c()
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    missing_packages <- c(missing_packages, pkg)
  }
}

# Report results
if (length(missing_packages) == 0) {
  cat("✓ All required packages are successfully installed!\n")
  cat("You can now run the main pipeline script.\n")
} else {
  cat("✗ The following packages are missing:\n")
  for (pkg in missing_packages) {
    cat("  -", pkg, "\n")
  }
  cat("Please install the missing packages manually.\n")
}

# =============================================================================
# 6. Load and Test Key Packages
# =============================================================================

cat("Testing package loading...\n")

tryCatch({
  library(TCGAbiolinks)
  library(DESeq2)
  library(org.Hs.eg.db)
  library(pheatmap)
  library(survminer)
  cat("✓ Core packages loaded successfully!\n")
}, error = function(e) {
  cat("✗ Error loading packages:", e$message, "\n")
})

cat("Installation script completed.\n")
cat("If all packages are installed successfully, you can run 'pipeline_documented.r'\n") 