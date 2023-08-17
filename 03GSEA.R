
# The code performs gene set enrichment analysis (GSEA) on gene expression data using MSigDB gene sets 
# and visualizes the enrichment results for specific pathways.

# 1. Load Libraries and Data:
#
# Load various libraries for performing enrichment analysis and visualization.
# Access MSigDB gene sets for different categories and subcategories.
# Load gene expression data or relevant information for analysis.
#
# 2. Preprocess Data:
#
# Convert gene symbols to ENTREZIDs using the org.Hs.eg.db database.
# Merge converted gene IDs with the original data, filtering out missing values.
# Prepare a list of gene fold changes and associate them with ENTREZIDs.
#
# 3. Gene Set Enrichment Analysis (GSEA):
#
# Perform Gene Set Enrichment Analysis (GSEA) using the gene fold changes and MSigDB gene sets.
# Calculate enrichment scores and results for each gene set.
#
# 4. Visualization and Analysis:
#
# Generate a CSV file to store the GSEA results.
# Visualize GSEA results using the gseaplot function for selected gene sets, displaying enrichment scores and p-values.
#
# 5. Batch Analysis and Visualization:
#
# Define lists of specific gene sets for further analysis.
# Plot GSEA results for the selected gene sets individually.
# 


library(clusterProfiler)#GO,KEGG,GSEA
library(GOplot)#弦图，弦表图，系统聚类图
library(DOSE)
library(ggnewscale)
library(topGO)#绘制通路网络图
library(circlize)#绘制富集分析圈图
library(ComplexHeatmap)#绘制图例
library(msigdbr)
library(ggrepel) # 标签相关
library(patchwork)
library(ggplot2)




# Access MSigDB gene sets for different categories and subcategories
db_H <- msigdbr(species = "Homo sapiens", category = "H") %>% 
  dplyr::select(gs_name, entrez_gene)

db_2K <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "KEGG") %>% 
  dplyr::select(gs_name, entrez_gene)

db_2B <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "BIOCARTA") %>% 
  dplyr::select(gs_name, entrez_gene)

db_5 <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "BP") %>% 
  dplyr::select(gs_name, entrez_gene)
# Combine all the obtained gene sets into a single data frame
db <- rbind(db_H,db_2K,db_2B,db_5)
head(db)



# Load gene expression data or other relevant data for GSEA analysis
info <- read.csv('./E746_A750del_VS_L858R.csv')

# Extract gene names from the loaded data
geneList <- info
genename <- as.character(geneList[,1])

# Convert gene symbols to ENTREZIDs using org.Hs.eg.db
library(org.Hs.eg.db)                        
geneList_tr <- bitr(genename,fromType = 'SYMBOL',toType = "ENTREZID",OrgDb = org.Hs.eg.db) 
colnames(geneList_tr)[1]<-"Gene"

# Merge the converted gene IDs with the original data
library('dplyr')
tb<-inner_join(geneList_tr,geneList,by = c("Gene"='X'))
tb<-na.omit(tb)
tb <- tb[order(-tb$log2FoldChange),]
geneList <- tb$log2FoldChange
names(geneList) <- as.character(tb$ENTREZID)
geneList

# Perform GSEA analysis using the obtained gene expression data and MSigDB gene sets
gsea <- GSEA(geneList,
             TERM2GENE=db,
             pvalueCutoff = 1
)
result <- gsea@result

# Perform gene set enrichment analysis using GO terms
Go_gseresult <- gseGO(geneList, 'org.Hs.eg.db', keyType = "ENTREZID", ont="all", nPerm = 1000, minGSSize = 10, maxGSSize = 1000, pvalueCutoff=1)

# Save the GSEA results to a CSV file
write.csv(gsea,'./gesa.csv')

# Plot GSEA results using gseaplot
gseaplot(Go_Reactomeresult,1,pvalue_table = TRUE)

# Define a list of specific gene sets for further analysis
list1 <- c("GOBP_DEVELOPMENTAL_CELL_GROWTH","GOBP_RIBOSOME_BIOGENESIS", "GOBP_HISTONE_LYSINE_DEMETHYLATION","GOBP_RIBONUCLEOPROTEIN_COMPLEX_BIOGENESIS") 

# Iterate through the list and create separate PDF plots for each gene set
for (path in list1){
  pdf(paste(path,".pdf"))
  a <- gseaplot(gsea,path, pvalue_table = TRUE,title = "") 
  print(a)
  dev.off()
}






