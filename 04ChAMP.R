#The code is a series of operations for analyzing DNA methylation data using the ChAMP package. 
#
# 1. Package Installation and Library Loading
#
# 2. Data Loading:
#   
# Mutation data (mut_df) and DNA methylation data (pd) are loaded.
# Row names are adjusted for the methylation data to match with the first column of pd.
# 
# 3. Data Preprocessing:
#   
# Unwanted columns are removed from the methylation data (pd).
# A grouping list (group_list) is extracted from the mutation data for future use.
#
# 4. Data Normalization:
#   
# Methylation data is normalized using champ.norm with specific parameters, and the normalized data is stored in myNorm.
# 
# 5. Differentially Methylated Positions (DMP) Analysis:
#   
# DMPs are identified using champ.DMP with specific parameters and are stored in myDMP.
# A subset of DMPs is extracted (Del19_to_p.L858R) and saved as a CSV file.
# 
# 6. DMP Visualization and Differentially Methylated Regions (DMR) Analysis:
#   
# The DMP.GUI function is used to visualize the identified DMPs based on the specified subset.
# DMRs are detected using champ.DMR based on the DMPs identified earlier (Del19_to_p.L858R).
# 
# 7. Gene Set Enrichment Analysis (GSEA):
#   
# Gene Set Enrichment Analysis is performed using champ.GSEA to identify enriched gene sets based on methylation data and other parameters.
# 
# 8. Gene Labeling and Sorting:
#   
# Genes from the DMP analysis results (result) are labeled as 'L858R' or '19Del' based on p-values and log fold changes.
# Genes are then sorted based on p-values and log fold changes.
#


# Set a timeout for BiocManager operations
options(timeout = 1200L)

# Install required Bioconductor packages
#BiocManager::install("ChAMPdata")
#BiocManager::install("ChAMP")

library(ChAMPdata)
library(ChAMP)


# Install required Bioconductor packages
mut_df <- read.csv("mut_df.csv")
pd <- read.csv("me_value.csv",header = T)
rownames(pd) <- pd$X
pd <- pd[,-1]
group_list <- mut_df$Amino_Acid_Change

# Normalize methylation data using ChAMP
myNorm <- champ.norm(beta=pd,arraytype="450K",cores=12) 

# Identify differentially methylated positions (DMPs) using ChAMP
myDMP <- champ.DMP(beta=pd,
                   pheno=group_list,
                   adjPVal = 1,
                   arraytype = "450K")
dim(myDMP$Del19_to_p.L858R)

# Extract a specific subset of DMPs and save the results to a CSV file
result <- myDMP$Del19_to_p.L858R
write.csv(result,"./Del19_to_p.L858R_raw.csv")

# Launch the ChAMP DMP GUI
DMP.GUI(myDMP$Del19_to_p.L858R,
        beta = myNorm,
        pheno=group_list)

# Identify differentially methylated regions (DMRs) using ChAMP
myDMR <- champ.DMR(myDMP$Del19_to_p.L858R,
                   beta = myNorm,
                   pheno=group_list,
                   adjPvalDmr=1)

# Perform Gene Set Enrichment Analysis (GSEA) using ChAMP
champ.GSEA(beta=myNorm,
           DMP=myDMP[[1]],
           CpGlist=NULL,
           Genelist=NULL,
           pheno=group_list,
           method="fisher",
           arraytype="450K",
           Rplot=TRUE,
           adjPval=1,
           cores=1) 

# Set the p-value threshold
pvalue <- 0.05

# Label significant genes based on log fold change and p-value
gene_diff <- myDMP$Del19_to_p.L858R
gene_diff$sig <- ifelse(gene_diff$P.Value < pvalue, 'L858R',
                        ifelse(gene_diff$P.Value < pvalue, '19Del', 'non-significant'))

# Sort the genes based on p-value and log fold change
gene_diff <- gene_diff[order(gene_diff$P.Value, gene_diff$logFC, decreasing = c(FALSE, TRUE)), ]

# Save differentially expressed genes (DEGs) to a CSV file
write.csv(as.data.frame(gene_diff), "DEGs_L858R_VS_p.Del19.csv")


# Select DEGs with specific labels for further analysis
gene_diff_select <- subset(gene_diff, sig %in% c('19Del', 'L858R'))
write.csv(gene_diff_select, file = "me_DEGs_select.csv", quote = FALSE)
L858R <- gene_diff %>% filter(gene_diff$sig=="L858R")
Del19 <- gene_diff %>% filter(gene_diff$sig=="19Del")
cpG_L858R <- CpG.GUI(CpG=rownames(L858R),arraytype="450K")
cpG_Del19 <- CpG.GUI(CpG=rownames(Del19),arraytype="450K")

# Quality control
champ.QC(beta = myNorm,
         pheno=group_list,
         mdsPlot=TRUE,
         densityPlot=TRUE,
         dendrogram=TRUE,
         PDFplot=TRUE,
         Rplot=TRUE,
         Feature.sel="None",
         resultsDir="./")

QC.GUI(beta=myNorm,
       pheno=group_list,
       arraytype="450K")





