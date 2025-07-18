# =============================================================================
# MRBP Breast Cancer Analysis Pipeline
# =============================================================================
# 
# This script performs comprehensive analysis of Mitochondrial RNA Binding 
# Proteins (MRBPs) in breast cancer using TCGA-BRCA transcriptome data.
#
# Analysis includes:
# - Data acquisition and preprocessing
# - Differential expression analysis
# - Survival analysis
# - Pathway enrichment analysis
# - Visualization (heatmaps, PCA, Venn diagrams)
#
# Author: Ritika Jaggi
# Date: May 2025
# =============================================================================

# =============================================================================
# SECTION 1: DATA ACQUISITION AND PREPARATION
# =============================================================================

# Load required libraries
library(TCGAbiolinks)
library(SummarizedExperiment)

# =============================================================================
# 1.1 Data Download from TCGA
# =============================================================================

# Query TCGA-BRCA transcriptome data
# This downloads gene expression quantification data from the TCGA breast cancer project
query <- GDCquery(
  project = "TCGA-BRCA",                    # TCGA Breast Cancer project
  data.category = "Transcriptome Profiling", # RNA-seq data
  data.type = "Gene Expression Quantification", # Gene expression counts
  workflow.type = "STAR - Counts"           # STAR alignment workflow
)

# Download the queried data
# Note: This will download a large amount of data (~ 6GB)
GDCdownload(query)

# Prepare the downloaded data into a SummarizedExperiment object
# This creates a structured object containing expression data and metadata
brca_data <- GDCprepare(query)

# =============================================================================
# 1.2 Data Processing and Filtering
# =============================================================================

# Extract expression counts and sample metadata
counts <- assay(brca_data)      # Gene expression count matrix
metadata <- colData(brca_data)  # Sample metadata (clinical information)

# Filter samples to include only primary tumor and normal tissue samples
# This ensures we're comparing cancer vs normal tissue
tumor_samples <- metadata[metadata$sample_type == "Primary Tumor",]
normal_samples <- metadata[metadata$sample_type == "Solid Tissue Normal",]

# Combine sample IDs for both groups
keep_samples <- c(rownames(tumor_samples), rownames(normal_samples))

# Subset the data to include only the filtered samples
counts <- counts[, keep_samples]
metadata <- metadata[keep_samples, ]

# =============================================================================
# SECTION 2: DIFFERENTIAL EXPRESSION ANALYSIS
# =============================================================================

library(DESeq2)

# =============================================================================
# 2.1 DESeq2 Analysis Setup
# =============================================================================

# Create DESeqDataSet object for differential expression analysis
# This object contains the count data and experimental design
dds <- DESeqDataSetFromMatrix(
  countData = counts,           # Gene expression count matrix
  colData = metadata,           # Sample metadata
  design = ~ sample_type        # Experimental design: compare by sample type
)

# Run DESeq2 analysis
# This performs normalization and statistical testing for differential expression
dds <- DESeq(dds)

# Extract results from the analysis
# Results are ordered by p-value (most significant first)
res <- results(dds)
res <- res[order(res$pvalue),]

# =============================================================================
# 2.2 Filtering Differentially Expressed Genes (DEGs)
# =============================================================================

# Filter for significantly differentially expressed genes
# Criteria: adjusted p-value < 0.05 AND |log2 fold change| > 1
degs <- res[which(res$padj < 0.05 & abs(res$log2FoldChange) > 1), ] 

# =============================================================================
# SECTION 3: MRBP GENE LIST FILTERING
# =============================================================================

# =============================================================================
# 3.1 Define MRBP Gene List
# =============================================================================

# List of Mitochondrial RNA Binding Proteins (MRBPs)
# These genes are involved in mitochondrial RNA metabolism, processing, and regulation
mrbp_list <- c("PUSL1","AURKAIP1","MRPL20","ATAD3B","ATAD3A","ACOT7","DNAJC11",
               "PARK7","SDHB","AKR7A2","AK2","MRPS15","YRDC","NDUFS5","OXCT2",
               "TRIT1","UQCRH","NSUN4","ATPAF1","NRDC","SCP2","CPT2","MRPL37",
               "ACOT11","PARS2","AK4","CRYZ","ACADM","KYAT3","ABCD3","DBT","AHCYL1",
               "ATP5PB","WARS2","MRPS21","TARS2","MRPL9","TDRKH","FLAD1","FDPS",
               "DAP3","SLC25A44","NAXE","MRPL24","PPOX","NDUFS2","SDHC","MGST3",
               "ADCY10","PRDX6","DARS2","MRPS14","TIMM17A","ANGEL2","LYPLAL1","IARS2",
               "MTARC1","PYCR2","COQ8A","MRPL55","GUK1","IBA57","COA6","TOMM20",
               "FH","TFB2M","RNASEH1","RDH14","HADHA","HADHB","MRPL33","TRMT61B",
               "LRPPRC","PREPL","MTIF2","PNPT1","NFU1","SPR","MTHFD2","MRPL53",
               "HTRA2","MRPL19","SUCLG1","PTCD3","IMMT","MRPL35","RPIA","MRPS5",
               "FAHD2A","COX5B","MRPL30","MRPS9","DBI","MMADHC","GPD2","FASTKD1",
               "METTL5","METTL8","METAP1D","MTX2","HIBCH","GLS","HSPD1","HSPE1",
               "MARS2","MAIP1","NIF3L1","NDUFB3","FASTKD2","CPS1","MRPL44","EFHD1",
               "NDUFA10","MTERF4","BOK","DTYMK","TRNT1","RPUSD3","CHCHD4","MRPS25",
               "OXNAD1","ACAA1","EXOG","HIGD1A","LARS2","DHX30","NME6","UQCRC1",
               "GPX1","HEMK1","CHDH","PDE12","PDHB","SUCLG2","NSUN3","CPOX","NIT2",
               "TOMM70","TRMT10C","ABHD10","TIMMDC1","NDUFB4","POLQ","ACAD9","MRPL3",
               "MRPS22","GFM1","MRPL47","DNAJC19","MCCC1","OPA1","ATP5ME","LETM1",
               "GRPEL1","LAP3","GUF1","OCIAD1","PAICS","NOA1","GRSF1","MTHFD2L","MRPL1","MRPS18C","NUDT9","PPA2","HADH","MGARP","MMAA","GATB","CASP3","SDHA","MRPL36","NDUFS6","NSUN2","FASTKD3","GOLPH3","OXCT1","NNT","MRPS30","SETD9","NDUFAF2","NLN","MRPS27","PTCD2","GFM2","DMGDH","COX7C","HSD17B4","ALDH7A1","HINT1","UQCRQ","VDAC1","SLC25A48","HSPA9","NDUFA2","HARS2","DELE1","MRPL22","THG1L","SFXN1","BPHL","FARS2","ALDH5A1","BAK1","TOMM6","MRPS10","MRPL2","MRPS18A","MRPL14","AARS2","MTO1","COX7A2","SMIM8","RARS2","LYRM2","RTN4IP1","QRSL1","MTRES1","HINT3","ECHDC1","MTHFD1L","MTRF1L","TFB1M","SOD2","MRPL18","MPC1","NDUFA4","MALSU1","CYCS","HIBADH","GARS1","MRPL32","MRPS24","OGDH","TBRG4","NIPSNAP2","MRPS17","CHCHD2","RCC1L","MDH2","GTPBP10","MTERF1","SLC25A13","PTCD1","ATP5MF-PTCD1","ATP5MF","PMPCB","DLD","PNPLA8","NDUFA5","SND1","CHCHD3","FMC1","ADCK2","MRPS33","GSTK1","FASTK","MICU3","GSR","PLPBP","VDAC3","LYPLA1","MRPL15","LACTB2","MRPS28","RMDN1","DECR1","PDP1","UQCRB","MTERF3","RIDA","COX6C","OXR1","MRPL13","NDUFB9","TOP1MT","CYC1","RECQL4","C8ORF82","AK3","STOML2","HINT2","GRHPR","ALDH1B1","RFK","AUH","MRPL50","NIPSNAP3A","HDHD3","STOM","NDUFA8","MRRF","PTRH1","PTGES2","TRUB2","ENDOG","MRPS2","PMPCA","MRPL41","IDI1","PITRM1","ATP5F1C","DHTKD1","NUDT5","YME1L1","MTPAP","TIMM23","OGDHL","CISD1","TFAM","DNA2","SUPV3L1","MCU","NUDT13","MRPS16","CHCHD1","VDAC2","PPIF","PRXL2A","GHITM","GLUD1","IDE","ALDH18A1","MRPL43","TWNK","SFXN3","PRDX3","OAT","ECHS1","MTG1","SIRT3","SLC25A22","MRPL17","HTATIP2","METTL15","CAT","PDHX","TIMM10","GLYAT","MRPL16","FTH1","PRDX5","ARL2","MRPL49","MRPL11","PC","NDUFV1","NDUFS8","CPT1A","MRPL21","CLPB","PDE2A","MRPL48","NDUFC2","NARS2","ACAT1","FDX1","DLAT","TIMM8B","REXO2","ATP5MG","RPUSD4","FOXRED1","NDUFA9","MRPL51","PHB2","MGST1","LDHB","MRPS35","ETFBKMT","DNM1L","YARS2","MYG1","CS","ATP5F1B","SHMT2","TSFM","MRPL42","SLC25A3","TXNRD1","ALDH1L2","MTERF2","MMAB","ALDH2","SIRT4","GATC","DIABLO","PUS1","PGAM5","MRPL57","MICU2","MTIF3","MRPS31","MTRF1","VWA8","SUCLA2","PCCA","CARS2","APEX1","METTL17","OXA1L","MRPL52","BCL2L2","SDR39U1","PRORP","SLC25A21","DMAC2L","TRMT5","ARG2","EXD2","ALDH6A1","DLST","GSTZ1","ALKBH1","SLIRP","CKMT1A","SQOR","DUT","LDHAL6B","LACTB","MTFMT","CLPX","COX5A","ETFA","IDH3A","MRPL46","MRPS11","POLG","IDH2","NGRN","MRPL28","NME4","MCRIP2","RHOT2","MRPS34","NDUFB10","ECI1","TRAP1","UQCRC2","EARS2","NDUFAB1","TUFM","GOT2","DUS2","DDX28","CYB5B","KARS1","CMC2","GCSH","COX4I1","ACSF3","GLOD4","MRM3","SLC25A11","C1QBP","ACADVL","SCO1","ELAC2","TTC19","PLD6","DHRS7B","POLDIP2","ERAL1","TEFM","LIG3","FKBP10","ACLY","COASY","COA3","SLC25A39","DCAKD","MRPL10","PHB","MRPL27","ACSF2","AKAP1","MRPS23","SEPTIN4","PTRH2","TACO1","FDXR","MRPL58","ATP5PD","MRPS7","MRPL38","OXLD1","MRPL12","SLC25A10","PYCR1","DCXR","FASN","NDUFV2","AFG3L2","OSBPL1A","ATP5F1A","ACAA2","ME2","RBFA","POLRMT","GPX4","ATP5F1D","TIMM13","MRPL54","MICOS13","LONP1","CLPP","ALKBH7","PET100","TIMM44","NDUFA7","MRPL4","QTRT1","PRDX2","GCDH","GADD45GIP1","TRMT1","MRPL34","GTPBP3","FKBP8","NDUFA13","UQCRFS1","COX6B1","TIMM50","TOMM40","BCAT2","BAX","ETFB","IDH3B","MRPS26","FASTKD5","MAVS","ACSS1","CDK5RAP1","RAB5IF","TOMM34","ATP5F1E","MTG2","MRPL39","ATP5PF","SOD1","ATP5PO","MRPS6","NDUFV3","GATD3A","YBEY","HDHD5","SLC25A18","BCL2L13","BID","MRPL40","COMT","SNAP29","NIPSNAP1","TST","MPST","PICK1","TOMM22","ACO2","CYB5R3","MCAT","TRMU","SCO2","MT-CO1","MT-CO2","MT-ATP6","MT-CO3","MT-CYB","GTPBP6","SLC25A6","HCCS","PDHA1","PRDX4","ACOT9","APOO","TIMM17B","HSD17B10","ABCB7","COX7B","PABPC5","TRMT2B","TIMM8A","ARMCX3","SLC25A53","SLC25A43","SLC25A5","AIFM1","SLC25A14","ABCD1","IDH3G","ACACA","ALDH9A1","ATP5MD","ATP5MPL","C12ORF65","CCDC58","CHCHD10","DNAJA3","ECH1","FAHD1","MCCC2","MRM1","MRPL23","MRPL45","MRPS12","MRPS18B","MRPS36","MTCH2","NDUFA6","NDUFS1","NDUFS3","NEU4","PAM16","PCK2","PRKACA","PTPMT1",
               "RDH13","SARS2","SLC25A24","SPRYD4","SSBP1","SURF1","TOP3A","VARS2") 

# =============================================================================
# 3.2 Map MRBP Gene Symbols to Ensembl IDs
# =============================================================================

library(org.Hs.eg.db)

# Convert MRBP gene symbols to Ensembl IDs for matching with TCGA data
# This is necessary because TCGA data uses Ensembl IDs while our MRBP list uses gene symbols
mrbp_ensembl <- mapIds(org.Hs.eg.db,
                       keys = mrbp_list,
                       column = "ENSEMBL",
                       keytype = "SYMBOL",
                       multiVals = "first")

# Remove any genes that couldn't be mapped (NA values)
mrbp_ensembl <- na.omit(mrbp_ensembl)

# =============================================================================
# 3.3 Filter Differential Expression Results for MRBPs
# =============================================================================

# Remove version numbers from Ensembl IDs in differential expression results
# TCGA data includes version numbers (e.g., ENSG00000000003.15)
# We need to remove these for proper matching
res_ids <- gsub("\\..*", "", rownames(res))

# Filter differential expression results to include only MRBPs
# Criteria: 
# 1. Gene is in our MRBP list
# 2. Gene is significantly differentially expressed (padj < 0.05)
# 3. Gene has substantial fold change (|log2FoldChange| > 1)
deg_mrbps <- res[
  res_ids %in% mrbp_ensembl &
    !is.na(res$padj) & res$padj < 0.05 &
    !is.na(res$log2FoldChange) & abs(res$log2FoldChange) > 1,
]

# =============================================================================
# SECTION 4: PCA AND CLUSTERING ANALYSIS
# =============================================================================

# =============================================================================
# 4.1 Principal Component Analysis (PCA)
# =============================================================================

# Perform variance stabilizing transformation (VST) for PCA
# VST normalizes the data and stabilizes variance across the mean
vsd <- vst(dds, blind = FALSE)

# Create PCA plot colored by sample type (tumor vs normal)
plotPCA(vsd, intgroup = "sample_type")

# =============================================================================
# 4.2 Gene ID to Symbol Conversion
# =============================================================================

library(org.Hs.eg.db)
library(AnnotationDbi)

# Remove version numbers from Ensembl IDs for annotation
ensembl_ids <- gsub("\\..*", "", rownames(assay(vsd)))

# Map Ensembl IDs to gene symbols for better interpretability
gene_symbols_vsd <- mapIds(org.Hs.eg.db, 
                           keys = ensembl_ids, 
                           column = "SYMBOL", 
                           keytype = "ENSEMBL", 
                           multiVals = "first")

# =============================================================================
# SECTION 5: HEATMAP VISUALIZATION
# =============================================================================

# =============================================================================
# 5.1 Prepare Heatmap Data
# =============================================================================

# Display summary of differentially expressed MRBPs
cat("Number of differentially expressed MRBPs:", nrow(deg_mrbps), "\n")
cat("Dimensions of differential expression results:", dim(deg_mrbps), "\n")

# Find MRBPs that are differentially expressed in our dataset
# This identifies which MRBPs from our list show significant changes in breast cancer
matched_genes <- rownames(deg_mrbps)[rownames(deg_mrbps) %in% rownames(assay(vsd))]

# Save the list of matched genes for reference
write.csv(matched_genes, file = "matched_genes.csv", row.names = FALSE)

# Extract expression data for differentially expressed MRBPs
deg_mrbps_ensg <- deg_mrbps[rownames(deg_mrbps) %in% rownames(assay(vsd)), ]
heatmap_data <- assay(vsd)[rownames(deg_mrbps_ensg), ]

# =============================================================================
# 5.2 Create Heatmaps
# =============================================================================

# Basic heatmap of MRBP expression
pheatmap::pheatmap(heatmap_data, 
                   cluster_rows = TRUE, 
                   cluster_cols = TRUE)

# Enhanced heatmap with custom color palette
pheatmap::pheatmap(heatmap_data, 
                   cluster_rows = TRUE, 
                   cluster_cols = TRUE, 
                   color = colorRampPalette(c("blue", "white", "red"))(50))

# Create a subset for better visualization (first 10 genes and 10 samples)
subset_data <- heatmap_data[1:10, 1:10]

# Heatmap of subset with smaller font sizes
pheatmap::pheatmap(subset_data, 
                   cluster_rows = TRUE, 
                   cluster_cols = TRUE, 
                   fontsize_row = 6,  
                   fontsize_col = 6)

# Row-normalized heatmap (each gene normalized to zero mean and unit variance)
pheatmap::pheatmap(subset_data, 
                   cluster_rows = TRUE, 
                   cluster_cols = TRUE, 
                   scale = "row")

# =============================================================================
# SECTION 6: SURVIVAL ANALYSIS
# =============================================================================

library(survival)
library(survminer)

# =============================================================================
# 6.1 Prepare Clinical Data
# =============================================================================

# Extract clinical information from the dataset
clinical <- colData(brca_data)

# =============================================================================
# 6.2 Example Survival Analysis for MRPL20
# =============================================================================

# Example analysis for one MRBP gene (MRPL20)
gene <- "MRPL20"

# Check if gene exists in the dataset
if (gene %in% rownames(counts)) {
  expression <- as.numeric(counts[gene, ])
} else {
  cat("Gene not found in the dataset\n")
}

# Find Ensembl ID for MRPL20
ensembl_id <- mapIds(org.Hs.eg.db, 
                     keys = "MRPL20", 
                     column = "ENSEMBL", 
                     keytype = "SYMBOL", 
                     multiVals = "first")
print(ensembl_id)

# =============================================================================
# 6.3 Expression Analysis and Visualization
# =============================================================================

# Clean gene ID and extract expression data
gene <- "ENSG00000242485"  # MRPL20 Ensembl ID
cleaned_gene <- sub("\\.\\d+$", "", gene)
cleaned_counts <- rownames(counts)
cleaned_counts <- sub("\\.\\d+$", "", cleaned_counts)

# Extract expression data for the gene
expression <- as.numeric(counts[36391, ])  # Using row index for MRPL20

# Create expression data frame
gene_expression_data <- data.frame(expression)

# Create boxplot comparing expression between tumor and normal samples
sample_conditions <- metadata$sample_type 
boxplot(expression ~ sample_conditions, 
        main = "Expression of MRPL20", 
        ylab = "Expression Level", 
        col = c("lightblue", "lightgreen"))

# =============================================================================
# 6.4 Survival Analysis Setup
# =============================================================================

# Get sample IDs from expression data
sample_ids <- colnames(counts)

# Match clinical data with expression data
clinical_matched <- clinical[sample_ids, ]

# Verify data alignment
cat("Expression data length:", length(expression), "\n")
cat("Clinical data rows:", nrow(clinical_matched), "\n")
cat("Data aligned:", length(expression) == nrow(clinical_matched), "\n")

# =============================================================================
# 6.5 Perform Survival Analysis
# =============================================================================

# Create high/low expression groups based on median expression
group <- ifelse(expression > median(expression), "High", "Low")

# Create survival object
# days_to_death: time to event
# vital_status == "Dead": event indicator (1 = event occurred, 0 = censored)
surv_object <- Surv(as.numeric(clinical_matched$days_to_death), 
                    clinical_matched$vital_status == "Dead")

# Fit survival model
fit <- survfit(surv_object ~ group)

# Create survival plot
ggsurvplot(fit, 
           data = clinical_matched, 
           pval = TRUE,           # Show p-value
           risk.table = TRUE)     # Show risk table

# =============================================================================
# SECTION 7: PATHWAY AND ENRICHMENT ANALYSIS
# =============================================================================

# =============================================================================
# 7.1 Install and Load Required Packages
# =============================================================================

# Install BiocManager if not already installed
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Install clusterProfiler for pathway analysis
BiocManager::install("clusterProfiler")

library(clusterProfiler)
library(org.Hs.eg.db)

# =============================================================================
# 7.2 Gene Ontology (GO) Enrichment Analysis
# =============================================================================

# Get Ensembl IDs of differentially expressed MRBPs
deg_ensembl <- gsub("\\..*", "", rownames(deg_mrbps))

# Convert Ensembl IDs to Entrez IDs for GO analysis
entrez_ids <- mapIds(
  org.Hs.eg.db,
  keys = deg_ensembl,
  column = "ENTREZID",
  keytype = "ENSEMBL",
  multiVals = "first"
)

# Perform GO enrichment analysis for all differentially expressed MRBPs
go_results_all_deg_mrbps <- enrichGO(
  gene          = entrez_ids,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",           # Biological Process
  pAdjustMethod = "BH",           # Benjamini-Hochberg correction
  pvalueCutoff  = 0.05,          # P-value threshold
  qvalueCutoff  = 0.2            # FDR threshold
)

# Create dotplot of GO enrichment results
dotplot(go_results_all_deg_mrbps)

# =============================================================================
# 7.3 GO Analysis for 66 Matched MRBPs
# =============================================================================

# Clean Ensembl IDs for the 66 matched genes
matched_genes_clean <- gsub("\\..*", "", matched_genes)

# Convert to Entrez IDs
entrez_ids <- mapIds(
  org.Hs.eg.db,
  keys = matched_genes_clean,
  column = "ENTREZID",
  keytype = "ENSEMBL",
  multiVals = "first"
)
entrez_ids <- na.omit(entrez_ids)

# Perform GO enrichment for the 66 matched MRBPs
go_results <- enrichGO(
  gene          = entrez_ids,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2
)

# Visualize results
dotplot(go_results)

# =============================================================================
# 7.4 Reactome Pathway Analysis
# =============================================================================

# Install ReactomePA if not already installed
if (!requireNamespace("ReactomePA", quietly = TRUE))
  BiocManager::install("ReactomePA")
library(ReactomePA)

# Perform Reactome pathway enrichment analysis
reactome_results <- enrichPathway(
  gene          = entrez_ids,
  organism      = "human",       # Human data
  pvalueCutoff  = 0.05,         # P-value threshold
  pAdjustMethod = "BH",          # Adjustment method
  qvalueCutoff  = 0.2,          # FDR threshold
  readable      = TRUE           # Convert ENTREZ back to gene symbols
)

# Display results
head(reactome_results)

# Create visualizations
dotplot(reactome_results, showCategory = 20)
cnetplot(reactome_results, categorySize = "pvalue", foldChange = NULL)

# =============================================================================
# SECTION 8: SAVE RESULTS
# =============================================================================

# Save all analysis results to CSV files
write.csv(as.data.frame(deg_mrbps), "DEG_MRBPs.csv")
write.csv(as.data.frame(go_results_all_deg_mrbps), "GO_Enrichment_MRBPs.csv")
write.csv(as.data.frame(go_results), 'GO_66.csv')
write.csv(as.data.frame(reactome_results), 'reactomePA.csv')

# =============================================================================
# SECTION 9: ADDITIONAL ANALYSIS - MRBP MATCHING
# =============================================================================

# =============================================================================
# 9.1 Match MRBPs to BRCA Data and Extract Results
# =============================================================================

library(org.Hs.eg.db)
library(DESeq2)
library(dplyr)

# Map MRBP gene symbols to Ensembl IDs
mrpb_ensembl <- mapIds(org.Hs.eg.db,
                       keys = mrbp_list,
                       column = "ENSEMBL",
                       keytype = "SYMBOL",
                       multiVals = "first")

# Clean Ensembl IDs in VST data to remove version suffix
cleaned_ensembl_ids <- sub("\\.\\d+$", "", rownames(vsd))

# Match MRBPs with BRCA dataset
names(cleaned_ensembl_ids) <- rownames(vsd)
ensembl_to_rowname <- cleaned_ensembl_ids[cleaned_ensembl_ids %in% mrpb_ensembl]

# Get actual rownames that match MRBPs
matching_rows <- names(ensembl_to_rowname)

# Extract expression data for MRBPs
mrbp_expr_matrix <- assay(vsd)[matching_rows, ]

# Extract differential expression data for MRBPs
res_df <- as.data.frame(res)
res_df$cleaned_id <- sub("\\.\\d+$", "", rownames(res_df))
mrbp_deg_results <- res_df[res_df$cleaned_id %in% mrpb_ensembl, ]

# Create heatmap of MRBP expression
library(pheatmap)
pheatmap(mrbp_expr_matrix, 
         show_rownames = TRUE, 
         show_colnames = FALSE,
         cluster_rows = TRUE, 
         cluster_cols = TRUE,
         main = "Expression of MRBPs in BRCA")

# Save MRBP differential expression results
write.csv(mrbp_deg_results, "MRBP_DE_results.csv")

# =============================================================================
# SECTION 10: VISUALIZATION - VENN DIAGRAMS
# =============================================================================

# =============================================================================
# 10.1 Pairwise Venn Diagram
# =============================================================================

# Clear previous plots
grid.newpage()

# Load required package
library(VennDiagram)

# Clean Ensembl IDs from TCGA data
tcga_genes <- sub("\\.\\d+$", "", rownames(vsd))

# Map MRBP symbols to Ensembl IDs
mrbp_ensembl <- mapIds(org.Hs.eg.db,
                       keys = mrbp_list,
                       column = "ENSEMBL",
                       keytype = "SYMBOL",
                       multiVals = "first")

# Remove NAs
mrbp_ensembl <- na.omit(mrbp_ensembl)

# Create pairwise Venn diagram
venn.plot <- draw.pairwise.venn(
  area1 = length(mrbp_ensembl),
  area2 = length(tcga_genes),
  cross.area = length(intersect(mrbp_ensembl, tcga_genes)),
  category = c("MRBPs", "TCGA-BRCA Genes"),
  fill = c("skyblue", "lightgreen"),
  lty = "blank",
  cex = 2,
  cat.cex = 1.5,
  cat.pos = c(-20, 20),
  cat.dist = 0.05
)

grid.draw(venn.plot)

# =============================================================================
# 10.2 Triple Venn Diagram
# =============================================================================

library(grid)

# Clear the plotting area
grid.newpage()

# Define the three sets for triple Venn diagram
# mrbp_ensembl: our MRBP gene list
# tcga_genes: genes in TCGA-BRCA dataset
# deg_genes: differentially expressed genes

# Clean versions if needed
mrbp_ensembl <- sub("\\.\\d+$", "", mrbp_ensembl)
tcga_genes <- sub("\\.\\d+$", "", rownames(counts))
deg_genes <- sub("\\.\\d+$", "", rownames(degs))

# Create triple Venn diagram
venn.plot <- draw.triple.venn(
  area1 = length(mrbp_ensembl),
  area2 = length(tcga_genes),
  area3 = length(deg_genes),
  n12 = length(intersect(mrbp_ensembl, tcga_genes)),
  n23 = length(intersect(tcga_genes, deg_genes)),
  n13 = length(intersect(mrbp_ensembl, deg_genes)),
  n123 = length(Reduce(intersect, list(mrbp_ensembl, tcga_genes, deg_genes))),
  category = c("MRBPs", "TCGA-BRCA", "DEGs"),
  fill = c("skyblue", "lightgreen", "plum"),
  lty = "blank",
  cex = 1.5,
  cat.cex = 1.3,
  cat.pos = c(-20, 20, 180),
  cat.dist = c(0.05, 0.05, 0.1)
)

grid.draw(venn.plot)

# =============================================================================
# 10.3 MRBP vs DEG Venn Diagram
# =============================================================================

# Map MRBP gene symbols to Ensembl IDs
ensembl_mrbps <- mapIds(org.Hs.eg.db, 
                        keys = mrbp_list, 
                        column = "ENSEMBL", 
                        keytype = "SYMBOL", 
                        multiVals = "first")

# Check mapping success
head(ensembl_mrbps)

# Clean Ensembl IDs from DEG results
deg_mrbps_cleaned <- sub("\\.\\d+$", "", rownames(deg_mrbps))

# Find overlapping genes
overlapping_genes <- intersect(ensembl_mrbps, deg_mrbps_cleaned)
print(overlapping_genes)

# Install and load VennDiagram if not already installed
if (!require(VennDiagram)) {
  install.packages("VennDiagram")
}
library(VennDiagram)

# Define the sets
mrbps_set <- ensembl_mrbps
degs_set <- rownames(deg_mrbps)

# Remove NA values
mrbps_set_clean <- na.omit(ensembl_mrbps)
degs_set_clean <- na.omit(rownames(deg_mrbps))

# Check data
cat("MRBPs (clean):", length(mrbps_set_clean), "\n")
cat("DEGs (clean):", length(degs_set_clean), "\n")

# Create Venn diagram
venn.plot <- venn.diagram(
  x = list(MRBPs = mrbps_set_clean, DEGs = degs_set_clean),
  category.names = c("MRBPs", "DEGs"),
  filename = NULL,
  output = TRUE,
  col = "black",
  fill = c("lightblue", "lightgreen"),
  alpha = 0.5,
  cex = 1.5,
  fontface = "bold",
  cat.cex = 1.5,
  cat.fontface = "bold",
  cat.col = c("blue", "green")
)

# Display the Venn diagram
grid.draw(venn.plot)

# =============================================================================
# SECTION 11: GENE NAME CONVERSION
# =============================================================================

# =============================================================================
# 11.1 Convert Matched Genes to Gene Names
# =============================================================================

# Install biomaRt if not already installed
if (!requireNamespace("biomaRt", quietly = TRUE)) {
  install.packages("BiocManager")
  BiocManager::install("biomaRt")
}

library(biomaRt)

# Remove version numbers from Ensembl gene IDs
ensembl_ids <- sub("\\..*", "", matched_genes)

# Connect to Ensembl database
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Query for gene symbols
genes_info <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters = "ensembl_gene_id",
  values = ensembl_ids,
  mart = mart
)

# Match back to original list
result <- merge(data.frame(ensembl_gene_id = ensembl_ids, original = matched_genes),
                genes_info, by = "ensembl_gene_id", all.x = TRUE)

# Display and save results
print(result)
write.csv(result, 'matched_gene_names.csv')

# =============================================================================
# END OF PIPELINE
# =============================================================================

cat("Pipeline completed successfully!\n")
cat("Output files have been saved to the working directory.\n")
cat("Key findings:\n")
cat("- Number of differentially expressed MRBPs:", nrow(deg_mrbps), "\n")
cat("- Number of matched genes:", length(matched_genes), "\n")
cat("- Survival analysis completed for MRPL20\n")
cat("- Pathway enrichment analysis completed\n")
cat("- All visualizations generated\n") 
