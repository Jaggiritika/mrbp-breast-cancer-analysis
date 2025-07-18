#data acquisition and preparation ####

library(TCGAbiolinks)

#fetching data
query <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)

#downloading data
GDCdownload(query)

#preparing data
brca_data <- GDCprepare(query)

#data processing ####

library(SummarizedExperiment)

#extract counts and metadata
counts <- assay(brca_data)
metadata <- colData(brca_data)

#filter to keep only primary tumor and normal tissue samples 
tumor_samples <- metadata[metadata$sample_type == "Primary Tumor",]
normal_samples <- metadata[metadata$sample_type == "Solid Tissue Normal",]
keep_samples <- c(rownames(tumor_samples), rownames(normal_samples))
counts <- counts[, keep_samples]
metadata <- metadata[keep_samples, ]

#differential expression analysis ####
library(DESeq2)

#create DESeqDataSet and run differential expression analysis
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = metadata,
                              design = ~ sample_type)
dds <- DESeq(dds)
res <- results(dds)
res <- res[order(res$pvalue),]

#filter DEGs based on adj. pvalue and log2fold change thresholds
degs <- res[which(res$padj < 0.05 & abs(res$log2FoldChange) > 1), ] 

#mrbp list filtering ####

#specify the mrbp genes list
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
               "GRPEL1","LAP3","GUF1","OCIAD1","PAICS","NOA1","GRSF1","MTHFD2L","MRPL1","MRPS18C","NUDT9","PPA2","HADH","MGARP","MMAA","GATB","CASP3","SDHA","MRPL36","NDUFS6","NSUN2","FASTKD3","GOLPH3","OXCT1","NNT","MRPS30","SETD9","NDUFAF2","NLN","MRPS27","PTCD2","GFM2","DMGDH","COX7C","HSD17B4","ALDH7A1","HINT1","UQCRQ","VDAC1","SLC25A48","HSPA9","NDUFA2","HARS2","DELE1","MRPL22","THG1L","SFXN1","BPHL","FARS2","ALDH5A1","BAK1","TOMM6","MRPS10","MRPL2","MRPS18A","MRPL14","AARS2","MTO1","COX7A2","SMIM8","RARS2","LYRM2","RTN4IP1","QRSL1","MTRES1","HINT3","ECHDC1","MTHFD1L","MTRF1L","TFB1M","SOD2","MRPL18","MPC1","NDUFA4","MALSU1","CYCS","HIBADH","GARS1","MRPL32","MRPS24","OGDH","TBRG4","NIPSNAP2","MRPS17","CHCHD2","RCC1L","MDH2","GTPBP10","MTERF1","SLC25A13","PTCD1","ATP5MF-PTCD1","ATP5MF","PMPCB","DLD","PNPLA8","NDUFA5","SND1","CHCHD3","FMC1","ADCK2","MRPS33","GSTK1","FASTK","MICU3","GSR","PLPBP","VDAC3","LYPLA1","MRPL15","LACTB2","MRPS28","RMDN1","DECR1","PDP1","UQCRB","MTERF3","RIDA","COX6C","OXR1","MRPL13","NDUFB9","TOP1MT","CYC1","RECQL4","C8ORF82","AK3","STOML2","HINT2","GRHPR","ALDH1B1","RFK","AUH","MRPL50","NIPSNAP3A","HDHD3","STOM","NDUFA8","MRRF","PTRH1","PTGES2","TRUB2","ENDOG","MRPS2","PMPCA","MRPL41","IDI1","PITRM1","ATP5F1C","DHTKD1","NUDT5","YME1L1","MTPAP","TIMM23","OGDHL","CISD1","TFAM","DNA2","SUPV3L1","MCU","NUDT13","MRPS16","CHCHD1","VDAC2","PPIF","PRXL2A","GHITM","GLUD1","IDE","ALDH18A1","MRPL43","TWNK","SFXN3","PRDX3","OAT","ECHS1","MTG1","SIRT3","SLC25A22","MRPL17","HTATIP2","METTL15","CAT","PDHX","TIMM10","GLYAT","MRPL16","FTH1","PRDX5","ARL2","MRPL49","MRPL11","PC","NDUFV1","NDUFS8","CPT1A","MRPL21","CLPB","PDE2A","MRPL48","NDUFC2","NARS2","ACAT1","FDX1","DLAT","TIMM8B","REXO2","ATP5MG","RPUSD4","FOXRED1","NDUFA9","MRPL51","PHB2","MGST1","LDHB","MRPS35","ETFBKMT","DNM1L","YARS2","MYG1","CS","ATP5F1B","SHMT2","TSFM","MRPL42","SLC25A3","TXNRD1","ALDH1L2","MTERF2","MMAB","ALDH2","SIRT4","GATC","DIABLO","PUS1","PGAM5","MRPL57","MICU2","MTIF3","MRPS31","MTRF1","VWA8","SUCLA2","PCCA","CARS2","APEX1","METTL17","OXA1L","MRPL52","BCL2L2","SDR39U1","PRORP","SLC25A21","DMAC2L","TRMT5","ARG2","EXD2","ALDH6A1","DLST","GSTZ1","ALKBH1","SLIRP","CKMT1A","SQOR","DUT","LDHAL6B","LACTB","MTFMT","CLPX","COX5A","ETFA","IDH3A","MRPL46","MRPS11","POLG","IDH2","NGRN","MRPL28","NME4","MCRIP2","RHOT2","MRPS34","NDUFB10","ECI1","TRAP1","UQCRC2","EARS2","NDUFAB1","TUFM","GOT2","DUS2","DDX28","CYB5B","KARS1","CMC2","GCSH","COX4I1","ACSF3","GLOD4","MRM3","SLC25A11","C1QBP","ACADVL","SCO1","ELAC2","TTC19","PLD6","DHRS7B","POLDIP2","ERAL1","TEFM","LIG3","FKBP10","ACLY","COASY","COA3","SLC25A39","DCAKD","MRPL10","PHB","MRPL27","ACSF2","AKAP1","MRPS23","SEPTIN4","PTRH2","TACO1","FDXR","MRPL58","ATP5PD","MRPS7","MRPL38","OXLD1","MRPL12","SLC25A10","PYCR1","DCXR","FASN","NDUFV2","AFG3L2","OSBPL1A","ATP5F1A","ACAA2","ME2","RBFA","POLRMT","GPX4","ATP5F1D","TIMM13","MRPL54","MICOS13","LONP1","CLPP","ALKBH7","PET100","TIMM44","NDUFA7","MRPL4","QTRT1","PRDX2","GCDH","GADD45GIP1","TRMT1","MRPL34","GTPBP3","FKBP8","NDUFA13","UQCRFS1","COX6B1","TIMM50","TOMM40","BCAT2","BAX","ETFB","IDH3B","MRPS26","FASTKD5","MAVS","ACSS1","CDK5RAP1","RAB5IF","TOMM34","ATP5F1E","MTG2","MRPL39","ATP5PF","SOD1","ATP5PO","MRPS6","NDUFV3","GATD3A","YBEY","HDHD5","SLC25A18","BCL2L13","BID","MRPL40","COMT","SNAP29","NIPSNAP1","TST","MPST","PICK1","TOMM22","ACO2","CYB5R3","MCAT","TRMU","SCO2","MT-CO1","MT-CO2","MT-ATP6","MT-CO3","MT-CYB","GTPBP6","SLC25A6","HCCS","PDHA1","PRDX4","ACOT9","APOO","TIMM17B","HSD17B10","ABCB7","COX7B","PABPC5","TRMT2B","TIMM8A","ARMCX3","SLC25A53","SLC25A43","SLC25A5","AIFM1","SLC25A14","ABCD1","IDH3G","ACACA","ALDH9A1","ATP5MD","ATP5MPL","C12ORF65","CCDC58","CHCHD10","DNAJA3","ECH1","FAHD1","MCCC2","MRM1","MRPL23","MRPL45","MRPS12","MRPS18B","MRPS36","MTCH2","NDUFA6","NDUFS1","NDUFS3","NEU4","PAM16","PCK2","PRKACA","PTPMT1","
              RDH13","SARS2","SLC25A24","SPRYD4","SSBP1","SURF1","TOP3A","VARS2") 

# Map MRBP gene symbols to Ensembl IDs
library(org.Hs.eg.db)
mrbp_ensembl <- mapIds(org.Hs.eg.db,
                       keys = mrbp_list,
                       column = "ENSEMBL",
                       keytype = "SYMBOL",
                       multiVals = "first")
mrbp_ensembl <- na.omit(mrbp_ensembl)
# Remove version numbers from Ensembl IDs in res
res_ids <- gsub("\\..*", "", rownames(res))
# Subset using Ensembl IDs
#deg_mrbps <- res[res_ids %in% mrbp_ensembl & !is.na(res$padj) < 0.05 & !is.na(abs(res$log2FoldChange) > 1, ]
deg_mrbps <- res[
  res_ids %in% mrbp_ensembl &
    !is.na(res$padj) & res$padj < 0.05 &
    !is.na(res$log2FoldChange) & abs(res$log2FoldChange) > 1,
]


#pca and clustering ####

#perform variance stabilizing transformation (VST) and plot pca
vsd <- vst(dds, blind = FALSE)
plotPCA(vsd, intgroup = "sample_type")

# Hierarchical clustering check DID nOT WORK ####
# debugging ####
#length(rownames(deg_mrbps)) 
#deg_mrbps <- degs[rownames(degs) %in% mrbp_list, ]


#converting ensembl ids to gene symbols ####
#BiocManager::install("org.Hs.eg.db")
library(org.Hs.eg.db)
library(AnnotationDbi)

#remove the version numbers(ENSG00000000003.15 -> ENSG00000000003)
ensembl_ids <- gsub("\\..*", "", rownames(assay(vsd)))

#mapping 
gene_symbols_vsd <- mapIds(org.Hs.eg.db, 
                           keys = ensembl_ids, 
                           column = "SYMBOL", 
                           keytype = "ENSEMBL", 
                           multiVals = "first")



#heatmap data selection and visualization ####

head(deg_mrbps)
dim(deg_mrbps)
# Get the gene symbols for rows in 'vsd'
#gene_symbols_vsd <- mapIds(org.Hs.eg.db, keys = rownames(assay(vsd)), column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")

# Subset the expression data (vsd) for the MRBP genes
#heatmap_data <- assay(vsd)[gene_symbols_vsd %in% deg_mrbps$symbol, ]
#head(rownames(assay(vsd)))
#keytypes(org.Hs.eg.db)



#table(is.na(gene_symbols_vsd))  # Check for any NA values (no mapping found)
#table(duplicated(gene_symbols_vsd))  # Check for duplicates (multiple gene symbols per Ensembl ID)

#gene_symbols_vsd <- gene_symbols_vsd[!is.na(gene_symbols_vsd)]
#gene_symbols_vsd <- gene_symbols_vsd[!duplicated(gene_symbols_vsd)]
#length(gene_symbols_vsd)

#ensure that the MRBP genes exist in the gene symbol list
#common_genes <- intersect(rownames(assay(vsd)), gene_symbols_vsd)
#length(common_genes)

#finding diff exp mrbps in vsd dataset (entire dataset with diff exp genes across tcga, ones which match >> gg)
head(rownames(assay(vsd)), 10)
deg_mrbps_ensg <- deg_mrbps[rownames(deg_mrbps) %in% rownames(assay(vsd)), ]
heatmap_data <- assay(vsd)[rownames(deg_mrbps_ensg), ]

matched_genes <- rownames(deg_mrbps)[rownames(deg_mrbps) %in% rownames(assay(vsd))] #66 found
write.csv(matched_genes, file = "matched_genes.csv", row.names = FALSE)


dim(heatmap_data)
pheatmap::pheatmap(heatmap_data, cluster_rows = TRUE, cluster_cols = TRUE)
pheatmap::pheatmap(heatmap_data, 
                   cluster_rows = TRUE, 
                   cluster_cols = TRUE, 
                   color = colorRampPalette(c("blue", "white", "red"))(50))


# Subset columns if needed (e.g., picking a group of interest)
subset_data <- heatmap_data[1:10,1:10]#For example, use the first 10 columns n rows aka first 10 tcga correlation with first 10 mrbps
pheatmap::pheatmap(subset_data, cluster_rows = TRUE, cluster_cols = TRUE)

pheatmap::pheatmap(subset_data, 
                   cluster_rows = TRUE, 
                   cluster_cols = TRUE, 
                   fontsize_row = 6,  # Smaller row labels
                   fontsize_col = 6)  # Smaller column labels



pheatmap::pheatmap(subset_data, 
                   cluster_rows = TRUE, 
                   cluster_cols = TRUE, 
                   scale = "row")  # Normalize each gene (row) to have zero mean and unit variance




####


#survival analysis ####
library(survival)
library(survminer) 

clinical <- colData(brca_data)  # check for vital_status, days_to_death, etc.

# Example for one MRBP gene
gene <- "MRPL20"

if (gene %in% rownames(counts)) {
  expression <- as.numeric(counts[gene, ])
} else {
  cat("Gene not found in the dataset\n")
}
head(rownames(counts), 10)
# Example for finding MRPL20 by Ensembl ID
ensembl_id <- mapIds(org.Hs.eg.db, keys = "MRPL20", column = "ENSEMBL", keytype = "SYMBOL", multiVals = "first")
print(ensembl_id)

gene <- "ENSG00000242485"
gene_in_counts <- gene %in% rownames(counts)
gene_in_counts

head(rownames(counts), 10)  # Print the first 10 row names
tail(rownames(counts), 10)  # Print the last 10 row names
cleaned_gene <- sub("\\.\\d+$", "", gene)  # Removes version number (e.g., .1, .2, etc.)
gene_in_counts <- cleaned_gene %in% rownames(counts)
gene_in_counts
cleaned_gene <- sub("\\.\\d+$", "", gene)
print(cleaned_gene)
cleaned_gene %in% rownames(counts)
head(rownames(counts), 10)

cleaned_counts <- rownames(counts)
cleaned_counts <- sub("\\.\\d+$", "", cleaned_counts)
cleaned_gene %in% cleaned_counts

expression <- as.numeric(counts[cleaned_gene, ])

head(rownames(counts), 10)  # check the first few row names of counts
head(cleaned_counts, 10)
cleaned_gene_versioned <- sub("\\.\\d+$", "", cleaned_gene)
cleaned_gene_versioned %in% rownames(counts)
grep(cleaned_gene_versioned, rownames(counts))



expression <- as.numeric(counts[36391, ])


head(expression)

gene_expression_data <- data.frame(expression)

# Make sure your metadata has the same number of rows as your counts data
sample_conditions <- metadata$sample_type 
boxplot(expression ~ sample_conditions, main = "Expression of MRPL20", ylab = "Expression Level")
boxplot(expression ~ sample_conditions, 
        main = "Expression of MRPL20", 
        ylab = "Expression Level", 
        col = c("lightblue", "lightgreen"))


length(expression)
nrow(clinical)
# Get sample IDs from expression (column names of counts matrix)
sample_ids <- colnames(counts)

# Subset clinical data to only those sample IDs
clinical_matched <- clinical[sample_ids, ]

# Now check dimensions again
length(expression) == nrow(clinical_matched)
#install.packages("xfun")
library("xfun")
#ggsurvplot(fit, data = clinical_matched, pval = TRUE, risk.table = TRUE)



# Assign groups based on expression levels
group <- ifelse(expression > median(expression), "High", "Low")

# Create survival object using the matched clinical data
surv_object <- Surv(as.numeric(clinical_matched$days_to_death), clinical_matched$vital_status == "Dead")

# Fit the model
fit <- survfit(surv_object ~ group)

# Plot the survival curves
ggsurvplot(fit, data = clinical_matched, pval = TRUE, risk.table = TRUE)


#pathway and enrichment analysis #####  

# Install BiocManager if not already installed
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Install clusterProfiler
BiocManager::install("clusterProfiler")

library(clusterProfiler)
library(org.Hs.eg.db)

deg_symbols <- rownames(deg_mrbps)
# Remove the version numbers
deg_ensembl <- gsub("\\..*", "", rownames(deg_mrbps))

entrez_ids <- mapIds(
  org.Hs.eg.db,
  keys = deg_ensembl,
  column = "ENTREZID",
  keytype = "ENSEMBL",
  multiVals = "first"
)
go_results_all_deg_mrbps <- enrichGO(
  gene          = entrez_ids,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2
)

dotplot(go_results_all_deg_mrbps)


#go for 66mrbps ####
# Remove everything after the "." in ENSEMBL IDs
matched_genes_clean <- gsub("\\..*", "", matched_genes)
entrez_ids <- mapIds(
  org.Hs.eg.db,
  keys = matched_genes_clean,
  column = "ENTREZID",
  keytype = "ENSEMBL",
  multiVals = "first"
)
entrez_ids <- na.omit(entrez_ids)


go_results <- enrichGO(
  gene          = entrez_ids,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2
)
dotplot(go_results)


# saving results ####
write.csv(as.data.frame(deg_mrbps), "DEG_MRBPs.csv")
write.csv(as.data.frame(go_results_all_deg_mrbps), "GO_Enrichment_MRBPs.csv")
write.csv(as.data.frame(go_results),'GO_66.csv')

#reactomePA ####
if (!requireNamespace("ReactomePA", quietly = TRUE))
  BiocManager::install("ReactomePA")
library(ReactomePA)

reactome_results <- enrichPathway(
  gene          = entrez_ids,
  organism      = "human",       # for human data
  pvalueCutoff  = 0.05,           # p-value threshold
  pAdjustMethod = "BH",           # adjust method
  qvalueCutoff  = 0.2,            # FDR threshold
  readable      = TRUE            # convert ENTREZ back to SYMBOL names
)
head(reactome_results)
dotplot(reactome_results, showCategory = 20)
cnetplot(reactome_results, categorySize = "pvalue", foldChange = NULL)
write.csv(as.data.frame(reactome_results),'reactomePA.csv')




















#Match MRBPs to BRCA data and extract results
# Step 1: Convert MRBP gene symbols to Ensembl IDs
library(org.Hs.eg.db)
library(DESeq2)
library(dplyr)

# Map SYMBOL to ENSEMBL
mrpb_ensembl <- mapIds(org.Hs.eg.db,
                       keys = mrbp_list,
                       column = "ENSEMBL",
                       keytype = "SYMBOL",
                       multiVals = "first")

# Step 2: Clean Ensembl IDs in counts/vsd to remove version suffix
cleaned_ensembl_ids <- sub("\\.\\d+$", "", rownames(vsd))

# Step 3: Match MRBPs with BRCA dataset
names(cleaned_ensembl_ids) <- rownames(vsd)  # temporarily label them
ensembl_to_rowname <- cleaned_ensembl_ids[cleaned_ensembl_ids %in% mrpb_ensembl]

# Get actual rownames that match MRBPs
matching_rows <- names(ensembl_to_rowname)

# Step 4: Extract expression data
mrbp_expr_matrix <- assay(vsd)[matching_rows, ]

# Step 5: Extract differential expression data
res_df <- as.data.frame(res)
res_df$cleaned_id <- sub("\\.\\d+$", "", rownames(res_df))
mrbp_deg_results <- res_df[res_df$cleaned_id %in% mrpb_ensembl, ]

# Optional: Heatmap of MRBP expression
library(pheatmap)
pheatmap(mrbp_expr_matrix, show_rownames = TRUE, show_colnames = FALSE,
         cluster_rows = TRUE, cluster_cols = TRUE,
         main = "Expression of MRBPs in BRCA")

# Optional: Save results
write.csv(mrbp_deg_results, "MRBP_DE_results.csv")



###

grid.newpage()  # Clears the previous plot (like heatmap)

# Load required package
library(VennDiagram)

# Clean Ensembl IDs from TCGA data (e.g., from vsd or counts rownames)
tcga_genes <- sub("\\.\\d+$", "", rownames(vsd))  # or rownames(counts) or rownames(res)

# Map MRBP symbols to Ensembl IDs
mrbp_ensembl <- mapIds(org.Hs.eg.db,
                       keys = mrbp_list,
                       column = "ENSEMBL",
                       keytype = "SYMBOL",
                       multiVals = "first")

# Remove NAs
mrbp_ensembl <- na.omit(mrbp_ensembl)

# Venn Diagram
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


library(grid)

# Clear the plotting area
grid.newpage()

# Assume these are all vectors of ENSEMBL IDs (same version format)
# mrbp_ensembl <- your list of MRBP Ensembl IDs
# tcga_genes <- rownames(counts) or rownames(vsd), depending on which you want
# deg_genes <- rownames(res_sig) from DESeq2 output

# Clean versions if needed
# mrbp_ensembl <- sub("\\.\\d+$", "", mrbp_ensembl)
# tcga_genes <- sub("\\.\\d+$", "", rownames(counts))
# deg_genes <- sub("\\.\\d+$", "", rownames(res_sig))

# Triple Venn
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



# Assuming your MRBP list is stored in `mrbps` (a vector of gene symbols)
library(org.Hs.eg.db)

# Map MRBP gene symbols to Ensembl IDs
ensembl_mrbps <- mapIds(org.Hs.eg.db, keys = mrbp_list, column = "ENSEMBL", keytype = "SYMBOL", multiVals = "first")
# Check if mapping was successful
head(ensembl_mrbps)
# Assuming DEGs are in `deg_mrbps` and they are in Ensembl format (rownames of the DEG dataset)
# Find overlap of MRBPs in DEGs

deg_mrbps_cleaned <- sub("\\.\\d+$", "", rownames(deg_mrbps))

# Now check for the overlap again
overlapping_genes <- intersect(ensembl_mrbps, deg_mrbps_cleaned)

# Print the overlapping genes
print(overlapping_genes)

# Install and load the necessary package if not already installed
if (!require(VennDiagram)) {
  install.packages("VennDiagram")
}
library(VennDiagram)

# Define the sets
mrbps_set <- ensembl_mrbps
degs_set <- rownames(deg_mrbps)
# Remove NA values from both MRBP and DEG datasets
mrbps_set_clean <- na.omit(ensembl_mrbps)
degs_set_clean <- na.omit(rownames(deg_mrbps))

# Check for NA values removal
length(mrbps_set_clean)  # Should be the length without NA
length(degs_set_clean)   # Should be the length without NA

# Create the Venn diagram
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



#converting matched_genes (66) to gene names ####
# Install biomaRt if not already installed
if (!requireNamespace("biomaRt", quietly = TRUE)) {
  install.packages("BiocManager")
  BiocManager::install("biomaRt")
}

library(biomaRt)

# Your vector of Ensembl gene IDs with version numbers
#matched_genes <- c("ENSG00000121691.7", "ENSG00000186642.16", "ENSG00000132837.15", "ENSG00000183010.17")

# Remove version numbers using sub
ensembl_ids <- sub("\\..*", "", matched_genes)

# Connect to Ensembl
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Query for gene symbols
genes_info <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters = "ensembl_gene_id",
  values = ensembl_ids,
  mart = mart
)

# Optional: match back to your original list
result <- merge(data.frame(ensembl_gene_id = ensembl_ids, original = matched_genes),
                genes_info, by = "ensembl_gene_id", all.x = TRUE)

# View result
print(result)
write.csv(result, 'matched_gene_names.csv')