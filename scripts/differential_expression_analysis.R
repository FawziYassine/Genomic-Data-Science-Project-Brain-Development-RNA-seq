# Differential Expression Analysis between the sample groups (adult vs. fetal)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install()

install.packages("dplyr")
install.packages("ggplot2")
install.packages("ggthemes")
install.packages("rlang")
remotes::install_github()
BiocManager::install("AnnotationDbi")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("DESeq")

library(SummarizedExperiment)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(DESeq2)
library(dplyr)
library(ggplot2)
library(ggthemes)
library(EnhancedVolcano)
setwd("~/brain")
par(pch = 19)

# read a merged_counts file
merged_counts = read.table("merged_counts.tsv", quote = "", sep = '\t')

# Read phenotype sample data
pheno_data = read.csv("phenotype_data.tsv", quote = "", sep = '\t')

# normalization: Read Per Million (RPM)  
x = as.matrix(merged_counts)
counts_RPM = t(t(x) * 1e6 / colSums(x))

table(pheno_data$sex)
table(pheno_data$age_group, useNA = "ifany")
table(pheno_data$sex, pheno_data$age_group)


# DESeq2 Differential Expression Analysis
# Create DESeq2 object 
deseq.dat = DESeqDataSetFromMatrix(countData = merged_counts, colData = pheno_data, design = ~ age_group + RIN)

# rownames(deseq.dat) = sub("\\.\\d+$", "", rownames(deseq.dat))

dds = DESeq(deseq.dat)

# DESeq2 results
res_deseq2 = results(dds, contrast = c("age_group", "adult", "fetal"))

head(res_deseq2)
res_deseq2_shrunk = lfcShrink(dds=dds, contrast = c("age_group", "adult", "fetal"), res = res_deseq2, type = "ashr")

# MA-plot: plot log2 fold-changes (on the y-axis) versus the mean of normalized counts (on the x-axis) for differentially expressed genes
pval_threshold = 10e-3
plotMA(res, alpha = pval_threshold)

# add gene symbol to res
gene_symbol = read.table("gencode.v28.symbols.txt", header = TRUE, na.strings = "n/a", col.names = c("gene", "symbol"))
table(duplicated(gene_symbol$symbol))
length(grep("RP4*", gene_symbol$symbol))
# check foxp2 gene
(foxp2.idx = which(gene_symbol$symbol == "FOXP2"))
res_deseq2_shrunk = as.data.frame(res_deseq2_shrunk)
res_deseq2_shrunk$row = rownames(res_deseq2_shrunk)
res_deseq2_shrunk_annotated = merge(res_deseq2_shrunk,gene_symbol, by.x = "row", by.y = "gene", all = T)
# check the number of annotated genes
sum(is.na(res_deseq2_shrunk_annotated$symbol))

# check foxp2 gene
foxp2.idx = which(res_deseq2_shrunk_annotated$symbol == "FOXP2")
res_deseq2_shrunk_annotated[foxp2.idx,]
# identify genes with FDR < 0.05
significant_differ_res_deseq2_shrunk = subset(res_deseq2_shrunk_annotated, res_deseq2_shrunk$padj < 0.05)
dim(significant_differ_res_deseq2_shrunk)
sum(is.na(significant_differ_res_deseq2_shrunk$symbol))

# Sort by increasing adjusted p-value
sorted_significant_differ_res_deseq2_shrunk = significant_differ_res_deseq2_shrunk[sort.list(significant_differ_res_deseq2_shrunk$padj,
                                                                                             decreasing = FALSE),]

# Volcano plot
EnhancedVolcano(sorted_significant_differ_res_deseq2_shrunk, FCcutoff = 1, pCutoff = pval_threshold,
                sorted_significant_differ_res_deseq2_shrunk$symbol,
                x = 'log2FoldChange', y = 'pvalue', 
                legend = c("Non-significant", "Passed log2 fold-change threshold", "Passed the p-value threshold", "Passed both thresholds"),
                legendPosition = 'right',
                legendLabSize = 9,
                legendIconSize = 2,
                widthConnectors = 0.2,
                colConnectors = 'grey30',
                colAlpha = 1)
res_deseq2_shrunk_annotated[which(res_deseq2_shrunk_annotated$symbol == "MT1B"),]


upload_df = data.frame(gene = sorted_significant_differ_res_deseq2_shrunk$row, name = sorted_significant_differ_res_deseq2_shrunk$symbol,
                       log2fc = sorted_significant_differ_res_deseq2_shrunk$log2FoldChange, pval = sorted_significant_differ_res_deseq2_shrunk$pvalue,
                       padj = sorted_significant_differ_res_deseq2_shrunk$padj)       
write.table(upload_df, "results.tsv", quote = FALSE, sep = '\t') 

# Principal Components Analysis (PCA) plot. A PCA plot shJows clusters of samples based on their similarity.
# It reduces the overwhelming number of dimensions by constructing principal components (PCs).
 
rld = rlogTransformation(dds)
plotPCA(rld, intgroup = "age_group") + ggtitle("PCA Plot of rLog Data") + theme_few()
