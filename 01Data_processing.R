# The code seems is for data processing, analysis, and conversion for genetic and epigenetic data. 
#   
# 1. Retrieve and Filter Data (expr and methylation matrix):
#   
# Gene-specific data is combined into a single dataframe using the retrieve_data function.
# Data is filtered based on specific criteria such as start and end thresholds, effects, and amino acid changes using the retrieve_and_filter_data function.
# Specific mutations like "Del19" and "L858R" are identified, and unique samples are selected.
# 
# 2. Data Integration and Transformation:
#   
# Methylation expression data (me) and gene expression data (expr) are loaded from files.
# Sample IDs are processed to match between data sets by modifying separators and finding common IDs.
# Data is filtered and matched using the sample IDs, and missing values are removed.
# 
# 3. ID Conversion:
#   
# Gene expression data (df) is converted from Ensembl IDs to gene IDs using a gene annotation file (gtf_data) for better interpretation.
# Methylation data (tb) is annotated with additional information from an annotation file (anno) based on unique IDs.
# Both converted gene expression data and annotated methylation data are saved into separate files.
#
# 4. Data Sorting and Export:
#   
# The mut_df dataframe is grouped and sorted by mutation type and sample IDs.
# Sample IDs are adjusted to match the sorted order.
# Gene expression and methylation data are exported to separate CSV files.


## Code begins:
## Set the working directory
# setwd("./data/02_exps/TCGA")


## Define a function to retrieve data from files for a given gene
retrieve_data <- function(file_paths, gene) {
  dbs <- list()
  
  for (i in 1:length(file_paths)) {
    file_path <- file_paths[i]
    db <- read.table(file_path, header = TRUE, sep = "\t")
    db <- db %>% filter(gene == gene)
    dbs[[i]] <- db
  }
  
  combined_db <- do.call(rbind, dbs)
  return(combined_db)
}

## Call retrieve_data to get gene-specific data
file_paths <- c(
  "./id_retrieve/TCGA-LUAD.varscan2_snv.tsv",
  "./id_retrieve/TCGA-LUAD.somaticsniper_snv.tsv",
  "./id_retrieve/TCGA-LUAD.mutect2_snv.tsv",
  "./id_retrieve/TCGA-LUAD.muse_snv.tsv"
)

gene_name <- "EGFR"
db <- retrieve_data(file_paths, gene_name)


## Define a function to retrieve and filter data based on certain criteria
## stop and end pos for 19Del sample_id
## L858R
## id:label
retrieve_and_filter_data <- function(db, start_threshold, end_threshold, effect_filter, amino_acid_change_filter) {
  Del19 <- db %>%
    filter(start > start_threshold & end < end_threshold) %>%
    filter(effect == effect_filter)
  Del19_uni <- subset(Del19, !duplicated(Sample_ID))
  
  L858R <- db %>%
    filter(Amino_Acid_Change == amino_acid_change_filter)
  L858R_uni <- subset(L858R, !duplicated(Sample_ID))
  
  mut_df <- rbind(Del19_uni, L858R_uni)
  mut_df <- mut_df[, c(1, 8)]
  rownames(mut_df) <- make.unique(mut_df$Sample_ID)
  
  return(mut_df)
}

## Call retrieve_and_filter_data function
start_threshold <- 55174701
end_threshold <- 55174840
effect_filter <- "inframe_deletion"
amino_acid_change_filter <- "p.L858R"

mut_df <- retrieve_and_filter_data(db, start_threshold, end_threshold, effect_filter, amino_acid_change_filter)

## Further filtering and labeling
T790M <- filter(db,Amino_Acid_Change=="p.T790M")
mut_df <- filter(mut_df,!(Sample_ID%in%intersect(mut_df$Sample_ID,T790M$Sample_ID)))
L858R <- db %>%
  filter(Amino_Acid_Change == "p.L858R")
L858R_uni <- subset(L858R, !duplicated(Sample_ID))
mut_df[1:23,2] <- "Del19"  # Adding a label


## Import methylation expression data matrix
me <- read.table("./TCGA-LUAD.methylation450.tsv",header=T,sep="\t")
row.names(me) <- me$Composite.Element.REF
me <- me[, -1]



## sample_id list
samp_id <- mut_df$Sample_ID
original_string <- samp_id 
samp_id <- gsub( "-","\\.", samp_id) #间隔符号转换
inter <- intersect(colnames(me),samp_id) #me表达矩阵和sample_id list共有的id

## Filter me rows using the ID list
me <- me[,inter]
me <- na.omit(me)
rownames(mut_df) <- samp_id
mut_df <- mut_df[inter,]
write.csv(me,"me_value.csv")
write.csv(mut_df,"mut_df.csv")

## Import gene expression matrix
expr <- read.table("./TCGA-LUAD.htseq_fpkm.tsv",header=T,sep="\t")
row.names(expr) <- expr$Ensembl_ID
expr <- expr[, -1]

## sample_id list
samp_id <- mut_df$Sample_ID
original_string <- samp_id 
samp_id <- gsub( "-","\\.", samp_id) #间隔符号转换
inter <- intersect(colnames(expr),samp_id) #expr表达矩阵和sample_id list共有的id

## Filter expr rows using the ID list
expr <- expr[,inter]
expr <- na.omit(expr)
rownames(mut_df) <- samp_id
mut_df <- mut_df[inter,]
write.csv(expr,"fpkm_ensembl.csv")
write.csv(mut_df,"mut_df.csv")


## =============================gene ENS id conversion=============================
## This section primarily involves data conversion and manipulation to match gene IDs with Ensembl IDs. 

## gtf_data: c("ensembl_id", "gene_id")
## Read the previously saved mut_df and df (gene expression data) files
mut_df <- read.csv("./data/02_exps/TCGA/mut_df.csv")
df <- read.csv("./data/02_exps/TCGA/fpkm_ensembl.csv")
rownames(df) <- df$X
df <- df[,-1]

## Specify the path to the gene annotation GTF file
gtf_path <- "./gencode.v43.chr_patch_hapl_scaff.annotation.gtf"

library("rtracklayer")
library("dplyr")

## Import data from the GTF file and extract necessary columns
gtf_data <- import(gtf_path)
gtf_data <- as.data.frame(gtf_data)
gtf_data <- data.frame(gtf_data$gene_id, gtf_data$gene_name)
colnames(gtf_data) <- c("ensembl_id", "gene_id")

## Remove duplicated gene IDs, keeping only the first occurrence
index <- duplicated(gtf_data[, 1])
gtf_data <- gtf_data[!index, ]
gtf_data$ensembl_id <- gsub("\\..*", "", gtf_data$ensembl_id)

## Check for duplicated entries in processed_gtf_data
table(duplicated(processed_gtf_data[, 1]))

## Prepare gene expression data (df) for conversion
cts <- round(df,2)
rownames(cts) <- gsub('\\..*','',rownames(cts))
cts$id <- rownames(cts)

## Convert Ensembl IDs to gene IDs
df <- left_join(cts,processed_gtf_data,by=c("id"="ensembl_id")) 
#table(duplicated(df$gene_id))
#head(df)
#index1 <- duplicated(df$gene_id)
df <- df[!duplicated(df$gene_id),]
#table(duplicated(df$gene_id))
df <- na.omit(df)
rownames(df) <- df$gene_id

library(dplyr)

## Sort sample IDs in mut_df according to the order in mut_df
mut_df <- mut_df %>%
  group_by(Amino_Acid_Change) %>%
  arrange(Amino_Acid_Change, Sample_ID)

samp_id <- mut_df$Sample_ID
samp_id <- gsub("-", "\\.", samp_id)

## Reorder df columns based on the sorted sample IDs
df <- df[, samp_id]

## Export the converted gene expression table
write.csv(df,'./geneid_fpkm.csv')
write.csv(mut_df,"mut_df.csv")


## =================================gene cg id convert==========================
## This section primarily involves methylation data annotations with corresponding gene information.


## Specify the path to the annotation file for DNA methylation data
anno_path <- "/./data/03_methylation/illuminaMethyl450_hg38_GDC.tsv"
tb <- read.csv("./data/03_methylation/TCGA/me_value.csv")
rownames(tb) <- tb$X
tb <- tb[,-1]

## Read the gene annotation file
anno <- read.table(anno_path)
colnames(anno) <- c("id","gene","chrom","chromStart","chromEnd","strand")
tb$id <- make.unique(rownames(tb))

## Merge DNA methylation data with annotation based on ID
tb_anno <- merge(tb,anno,by="id")
write.csv(tb_anno,"./me_anno.csv")

rownames(tb_anno) <- make.unique(tb_anno$gene)

## Remove unnecessary columns
me_anno <- me_anno[,c(-43,-42,-41,-40,-39,-1)]

## Save results
write.csv(tb_anno,"me_anno.csv")

## Sort sample IDs in mut_df according to the order in mut_df
mut_df <- mut_df %>%
  group_by(Amino_Acid_Change) %>%
  arrange(Amino_Acid_Change, Sample_ID)

samp_id <- mut_df$Sample_ID
samp_id <- gsub("-", "\\.", samp_id)

## Reorder me columns based on the sorted sample IDs
me <- me[, samp_id]

## Export the converted DNA methylation table and the mut_df table
write.csv(mut_df,"./me/mut_df.csv")
write.csv(me,"./me/me_value1.csv")

