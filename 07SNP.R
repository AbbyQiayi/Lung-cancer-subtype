# SNP analysis

# Load required packages
require(maftools)
# Set an option to treat strings as factors
options(stringsAsFactors = F) 
library(data.table)

# Load MAF data files using 'fread' function from the 'data.table' package
tmp1=fread('TCGA.LUAD.somaticsniper.maf')
tmp2=fread('TCGA.LUAD.muse.maf')
tmp3=fread('TCGA.LUAD.verscan.maf')
tmp4=fread('TCGA.LUAD.mutect.maf')

# Combine data frames using 'rbind'
df <- rbind(tmp4,tmp3,tmp2,tmp1)

# Write the combined data frame to a file named 'df.maf'
write.table(df,'df.maf',row.names = F,col.names = F)

# Read mutation data from a CSV file and remove the first column
mut_df <- read.csv("./data/02_exps/TCGA/mut_df.csv",header = T,sep = ',')
mut_df <- mut_df[,-1]

# Filter rows with specific amino acid changes
del19 <- mut_df %>% filter(mut_df$Amino_Acid_Change == 'Del19')
L858R <- mut_df %>% filter(mut_df$Amino_Acid_Change == 'p.L858R')


# Extract Sample IDs from filtered data
samp_19 <- del19$Sample_ID
samp_L858R <- L858R$Sample_ID

# Filter the main data frame based on extracted Sample IDs
db_19 <- df %>% filter(substr(df$Tumor_Sample_Barcode, 1, 16) %in% samp_19)
db_858 <- df %>% filter(substr(df$Tumor_Sample_Barcode, 1, 16) %in% samp_L858R)

## maftools
colnames(df)[colnames(df) == "B"] <- "NewColumn"

# Read clinical data from a TSV file using 'fread'
clinicalDat <- fread("TCGA-LUAD.GDC_phenotype.tsv")
all_tb <- read.maf(maf = df,clinicalData=clinicalDat,
               vc_nonSyn=names(tail(sort(table(db_19$Variant_Classification)))))

# Perform MAF analysis using 'read.maf' and store the results in 'obj' and 'obj_19'
obj_19 = read.maf(maf = db_19,clinicalData=clinicalDat,
                     vc_nonSyn=names(tail(sort(table(db_19$Variant_Classification)))))

obj_858 = read.maf(maf = db_858,clinicalData=clinicalDat,
                   vc_nonSyn=names(tail(sort(table(db_858$Variant_Classification)))))

tmp4_maf <- read.maf(maf = tmp4,clinicalData=clinicalDat,
                   vc_nonSyn=names(tail(sort(table(db_19$Variant_Classification)))))

# Create a PDF file for oncoplot and generate an oncoplot
pdf(file="oncoplot.pdf",width=10,height=10)
oncoplot(maf = obj, top = 10) # Plot the top 10 most frequent mutated genes
dev.off()

# Create a PDF file for plotmafSummary
getFields(obj)
pdf(file="plotmafSummary.pdf",width=10,height=10)
plotmafSummary(maf = obj, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
dev.off()


obj.titv = titv(maf = obj, plot = FALSE, useSyn = TRUE)
# Create a PDF file for Transition and Transversions summary
pdf(file="plotTiTv.pdf",width=8,height=6)
plotTiTv(res = obj.titv)
dev.off()

# Create a PDF file for somaticInteractions
pdf(file="somaticInteractions.pdf",width=10,height=10)
somaticInteractions(maf = obj, top = 20, pvalue = c(0.05, 0.1))
dev.off()

# Create a PDF file for forestPlot
DEL19 = read.maf(maf = db_19)
L858R = read.maf(maf = db_858)
#Considering only genes which are mutated in at-least in 5 samples in one of the cohort to avoid bias due to genes mutated in single sample.
L858R.vs.Del19 <- mafCompare(m1 = L858R, m2 = DEL19, m1Name = 'L858R', m2Name = 'Del19', minMut = 5)
print(L858R.vs.Del19)
#ForestPlot
pdf(file="forestPlot.pdf",width=9,height=6)
forestPlot(mafCompareRes = L858R.vs.Del19, pVal = 0.2)
dev.off()

# Create a PDF file for coOncoplot
genes = c("RBM10")
pdf(file="coOncoplot.pdf",width=10,height=4)
coOncoplot(m1 = L858R, m2 = DEL19, m1Name = 'L858R', m2Name = 'Del19', genes = genes, removeNonMutated = TRUE)
dev.off()

# Create a PDF file for coBarplot
pdf(file="coBarplot.pdf",width=6,height=6)
coBarplot(m1 = L858R, m2 = DEL19, m1Name = 'L858R', m2Name = 'Del19')
dev.off()

fab.ce = clinicalEnrichment(maf = all_tb)



