# 
# Read data
df <- read.csv("./me_value.csv", header = TRUE, sep = ",")
rownames(df) <- df$X
df <- df[, -1]
df <- na.omit(df)

# Read mutations
mut_df <- read.table("./mut_df.csv", header = TRUE, sep = ",")

# Prepare coldata
coldata <- data.frame(condition = mut_df$Amino_Acid_Change)

df <- me
# Run DESeq2 analysis
library("DESeq2")
dds <- DESeqDataSetFromMatrix(countData = round(df), colData = coldata, design = ~ condition)
dds <- DESeq(dds)
resultsNames(dds)
result <- results(dds)
result <- data.frame(result)[complete.cases(result),]

# Save results
write.csv(result, './fpkm_p.L858R_VS_Del19_deseq2.csv', row.names = TRUE)


# Set the p-value threshold
pvalue <- 0.05

# Filter and label significant genes
gene_diff <- result
gene_diff$sig <- ifelse(gene_diff$log2FoldChange > 0 & gene_diff$pvalue < pvalue, 'L858R',
                        ifelse(gene_diff$log2FoldChange < 0 & gene_diff$pvalue < pvalue, '19Del', 'non-significant'))

# Sort the genes
gene_diff <- gene_diff[order(gene_diff$pvalue, gene_diff$log2FoldChange, decreasing = c(FALSE, TRUE)), ]

# Save results to a CSV file
write.csv(as.data.frame(gene_diff), "DEGs_L858R_VS_p.Del19_FPKM_deseq2.csv")


#select DEGs
gene_diff_select <- subset(gene_diff, sig %in% c('19Del', 'L858R'))
write.csv(gene_diff_select, file = "fpkm_select_deseq2.csv", sep = '\t', col.names = NA, quote = FALSE)



# edgeR

library(edgeR)
library(limma)
library(statmod)

# Set the p-value threshold
pvalue <- 0.05

# Function to preprocess the data

design <- coldata$condition
dgelist <- DGEList(counts = df, group = design)

# Filter low count data using CPM normalization
keep <- rowSums(cpm(dgelist) > 1) >= 2
dgelist <- dgelist[keep, , keep.lib.sizes = FALSE]

# Normalize using TMM method
dgelist_norm <- calcNormFactors(dgelist, method = 'TMM')


# Function to perform differential expression analysis

design <- model.matrix(~design)

# Estimate dispersions
dge <- estimateDisp(dgelist_norm, design, robust = TRUE)

# Fit the model
fit <- glmFit(dge, design, robust = TRUE)

# Perform likelihood ratio test
lrt <- topTags(glmLRT(fit), n = nrow(dgelist_norm$counts))
lrt <- lrt$table

# Label significant genes
gene_diff <- lrt[order(lrt$PValue, lrt$logFC, decreasing = c(FALSE, TRUE)), ]
gene_diff$sig <- ifelse(gene_diff$log2FoldChange > 0 & gene_diff$pvalue < pvalue,  'L858R',
                        ifelse(gene_diff$log2FoldChange < 0 & gene_diff$pvalue < pvalue,'19Del', 'non-significant'))



# Save results to a CSV file
write.csv(as.data.frame(gene_diff), "p.L858R_VS_Del19_edgeR.csv")

# Output select
gene_diff_select <- subset(gene_diff, sig %in% c('19Del', 'L858R'))
write.csv(gene_diff_select, file = "select.csv", sep = '\t', col.names = NA, quote = FALSE)
result <- gene_diff


## volcano plot
result <- read.csv("/home/qiayi/R/shuoshi/data/Xena/me/me_DEGs_select.csv")

result$label <- ifelse(result$P.Value < 0.05 & abs(result$logFC) >0,rownames(result),"")
table(result$label)

library(ggplot2)
library(ggrepel)
ggplot(result, aes(x = log2FoldChange, y = -log10(pvalue), colour=sig)) +
  geom_point(alpha=0.4, size=3)+
  scale_color_manual(values=c("#ff4757","#546de5","#d2dae2"))+
  # 辅助线
  geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.4) +
  geom_hline(yintercept = -log10(cut_off_pvalue),
             lty=4,col="black",lwd=0.4) +
  # 坐标轴
  labs(x="log2(Fold Change)",
       y="-log10 (pvalue)")+
  theme_bw()+
  ggtitle("Volcano Plot")+
  # 图例
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position="right", 
        legend.title = element_blank()
  )+ 
  geom_text_repel(data = result, aes(x = logFC, y = -log10(PValue), label = label),size=2.7)
