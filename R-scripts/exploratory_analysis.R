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
devtools::install_github('alyssafrazee/RSkittleBrewer')

library(SummarizedExperiment)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(DESeq2)
library(dplyr)
library(ggplot2)
library(gplots)
library(ggthemes)
library(EnhancedVolcano)
setwd("~/brain")
par(pch = 19) 
tropical = c("darkorange", "dodgerblue", "hotpink", "limegreen", "yellow")
palette(tropical)

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

summary(counts_RPM)
sum(pheno_data$age_group == " ")
is.na(counts_RPM[1,])
sum(is.na(counts_RPM))

# Make the distribution of NA's by genes
gene_na = rowSums(is.na(counts_RPM))
gene_na[5]
gene_na[1:5]

# Make the distribution of NA's by samples
sample_na = colSums(is.na(counts_RPM))
sample_na[6]
table(sample_na)

# boxplot 1st sample with log2 transformation of the data
boxplot(log2(counts_RPM[,1]+1), col = 2, xlab = "sample",  ylab = "log2(counts[,1] + 1)")

# boxplot all samples without any transformation of the data
boxplot(counts_RPM, , col = 2, xlab = "samples", ylab = "counts")

# boxplot all samples with log2 transformation of the data
boxplot(log2(counts_RPM+1), col = 2, xlab = "samples",  ylab = "log2(counts + 1)")

# boxplot all samples with log2 transformation of the filtered data
counts_RPM = as.data.frame(counts_RPM)
fil_counts_RPM = filter(counts_RPM, rowMeans(counts_RPM) > 1)
dim(fil_counts_RPM)
boxplot(as.matrix(log2(fil_counts_RPM+1)), col=2, xlab = "samples",  ylab = "log2(filtered_counts + 1)")

# histogram: show the probability/frequency distribution of counts data
par(mfrow=c(1,2))
hist(log2(counts_RPM[,1]+1), col = 2)

# boxplot all samples with log2 transformation of the filtered data
hist(log2(fil_counts_RPM[,1]+1), col = 2, Xlab = "log2(filtered_counts[,1] + 1)")
low = log2(counts_RPM[,1]+1) < 1                                                               
head(low)
table(low)
hist(log2(counts_RPM[,2]+1), col = 2)
par(mfrow=c(1,1))
  
# densities
plot(density(log2(counts_RPM[,1]+1)), col = 2)
lines(density(counts_RPM[,2], col = 3))

# MA plot between 2 samples     
aa = log2(counts_RPM[,1]+1) + log2(counts_RPM[,2]+1)
mm = log2(counts_RPM[,1]+1) - log2(counts_RPM[,2]+1)
plot(aa, mm, col=2)
str(counts_RPM)

# Confirm that male samples have more genes on chromosome Y than females
rownames(counts_RPM) = sub("\\.\\d+$", "", rownames(counts_RPM))

chr = AnnotationDbi::select(org.Hs.eg.db, keys = rownames(counts_RPM), keytype = "ENSEMBL", columns = "CHR")
chr = chr[!duplicated(chr[,1]),]

# Confirm that the annotations still have the same sort as the counts
all(chr[,1] == rownames(counts_RPM))

# Select the chromosome Y samples
dim(chr[which(chr$CHR == "Y"),])
fil2_counts_RPM = filter(counts_RPM, chr$CHR == "Y")

# Male samples have more genes on chromosome Y than females
boxplot(colSums(fil2_counts_RPM) ~ pheno_data$sex)
points(colSums(fil2_counts_RPM) ~ jitter(as.numeric(pheno_data$sex)), col = as.numeric(pheno_data$sex))

# Create DESeq2 object 
deseq.dat = DESeqDataSetFromMatrix(countData = merged_counts, colData = pheno_data, design = ~ age_group + RIN)

# rownames(deseq.dat) = sub("\\.\\d+$", "", rownames(deseq.dat))
head(assay(deseq.dat))

dds = DESeq(deseq.dat)

# Principl Components Analysis (PCA) plot. A PCA plot shJows clusters of samples based on their similarity.
# It reduces the overwhelming number of dimensions by constructing principal components (PCs).
rld = rlogTransformation(dds)
plotPCA(rld, intgroup = "age_group") + ggtitle("PCA Plot of rLog Data") + theme_few()

# Heatmap: gives us an overview over similarities and dissimilarities between samples.
fil3_counts_RPM = filter(counts_RPM, rowMeans(counts_RPM) > 1)
dim(fil3_counts_RPM)
heatmap(as.matrix(fil3_counts_RPM))
coloramp =  colorRampPalette(c(3, "white", 2))(9)
heatmap(as.matrix(fil3_counts_RPM), col = coloramp, Rowv = NA, Colv = NA)
# heatmap of the expression data itself
heatmap.2(as.matrix(fil3_counts_RPM), col = coloramp, Rowv = NA, Colv = NA,
          dendogram="none", scale="row", trace="none")
#heatmap of the distances between samples
fil3_counts_RPM = log2(fil3_counts_RPM + 1)
dist.samples = dist(t(fil3_counts_RPM))
heatmap.2(as.matrix(dist.samples), col = coloramp, Rowv = NA, Colv = NA,
          dendogram="none", scale="row", trace="none")

# The rotation matrix contains the principal component loading. Each column of the rotation matrix contains the principal component loading vector.
# The component loading can be represented as the correlation of a particular variable on the respective PC(principal component).
pca_res2$rotation

