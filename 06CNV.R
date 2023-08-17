# CNV analysis
# ============ Preprocess data for GISTIC2.0 analysis ==========================

# Load required libraries
library(dplyr)         # For data manipulation
library(GISTIC)        # For GISTIC analysis
library(MutationAnnotationFormat)  # For handling mutation annotation files

# Read CNV data and mutation data
db <- read.table("./data/01_mut/CNV/luad_tcga_pub_segments.seg",header = T)
mut_df <- read.csv("./data/02_exps/TCGA/mut_df.csv",header = T,sep = ',')
mut_df <- mut_df[,-1]

# Filter mutation data for specific mutations
del19 <- mut_df %>% filter(mut_df$Amino_Acid_Change == 'Del19')
L858R <- mut_df %>% filter(mut_df$Amino_Acid_Change == 'p.L858R')
samp <- L858R$Sample_ID

# Filter CNV data based on selected samples
df <- db %>% filter(paste0(db$ID,'A') %in% samp)

# Write the filtered CNV data to a file
write.table(df,"L858R.seg",row.names = FALSE,col.names = FALSE,sep = '\t')

# Find the intersection of unique IDs and selected samples
intersect(paste0(unique(db$ID),'A'),samp)

## ============== gistic analysis of data processed by GISTIC2.0 ===============


# gistic analysis
lung.gistic <- readGistic(gisticAllLesionsFile="./data/01_mut/CNV/Del19/del19/del19/all_lesions.conf_75.txt",gisticAmpGenesFile="/home/qiayi/R/shuoshi/data/01_mut/CNV/Del19/del19/del19/amp_genes.conf_75.txt",gisticDelGenesFile="/home/qiayi/R/shuoshi/data/01_mut/CNV/Del19/del19/del19/del_genes.conf_75.txt",gisticScoresFile="/home/qiayi/R/shuoshi/data/01_mut/CNV/Del19/del19/del19/scores.gistic")
lung.gistic

# Write GISTIC summary to files
write.GisticSummary(gistic=lung.gistic,basename="plot/lung/L858R")

# Create a PDF file for GISTIC chromosomal plot
pdf(file="plot/gistic.pdf",width=10,height=6)
gisticChromPlot(gistic=lung.gistic,markBands="all")
dev.off()

# Create a PDF file for GISTIC bubble plot
pdf(file="plot/gisticBubble.pdf",width=4,height=4)
gisticBubblePlot(gistic=lung.gistic)
dev.off()

# Read mutation annotation file
lungmaf <- read.maf("./project/lung808/backup/03_mutation/test/total_sample/all_sample_mutations_filter.maf")
pdf(file="plot/cnv_oncoplot.pdf",width=4,height=4)
gisticOncoPlot(gistic=lung.gistic)
dev.off()

