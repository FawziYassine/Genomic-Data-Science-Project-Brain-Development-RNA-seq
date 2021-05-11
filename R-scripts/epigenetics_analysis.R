# Epigenetics and Differential Expression Analysis 

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("AnnotationHub")
BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene")
BiocManager::install("EnsDb.Hsapiens.v86")
BiocManager::install("mygene")
install.packages("httr")
packageurl = "https://cran.r-project.org/src/contrib/Archive/httr/httr_1.0.0.tar.gz"
install.packages(packageurl, repos = NULL, type = "source")

library(AnnotationHub)
library(AnnotationDbi)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
library(EnsDb.Hsapiens.v86)
library(mygene)
 
# load list of of differentially expressed gene
diff.expressed.genes = read.table("Differentially_Expressed_Genes-results.tsv", quote = "", sep = '\t')

# check to see if FXN was in these differentially expressed genes
diff.expressed.genes[which(diff.expressed.genes$name == "FXN"), ]

# AnnotationHub
ah = AnnotationHub()
ah
length(ah)
ah = subset(ah, species == "Homo sapiens")

# get the narrow peaks data (promoter associated histone modification H3K4me3) for the fetal brain cell line
fetal.brain = AnnotationHub::query(ah, c("EpigenomeRoadMap", "H3K4me3", "E081"))
fetal.brain.gr = fetal.brain[[2]]

# get the narrow peaks data (promoter associated histone modification H3K4me3) for the adult brain cell line
adult.brain = AnnotationHub::query(ah, c("EpigenomeRoadMap", "H3K4me3", "E073"))
adult.brain.gr = adult.brain[[2]]

# get the narrow peak data (promoter associated histone modification H3K4me3) for liver cell line
liver.line = AnnotationHub::query(ah, c("EpigenomeRoadMap", "H3K4me3", "Liver"))
liver.line.gr = liver.line[[2]]

# get the known genes from Tx database
txdb = TxDb.Hsapiens.UCSC.hg19.knownGene
txdb.genes = genes(txdb)

# check the str ucture of FXN gene (Entrez ID 2395) 
(fxn.gr = txdb.genes[which(txdb.genes$gene_id == "2395"), ])
subsetByOverlaps(genes(txdb), fxn.gr)
subsetByOverlaps(transcripts(txdb), fxn.gr)
subsetByOverlaps(exons(txdb), fxn.gr)
subsetByOverlaps(cds(txdb), fxn.gr)

# change "ENSEMBL" ID to "ENTREZID" that matches gene_id
diff.expressed.genes$gene = sub("\\.\\d+$", "", as.character(diff.expressed.genes$gene))
AnnotationDbi::keytypes(EnsDb.Hsapiens.v86)
# diff.expressed.genes = diff.expressed.genes[!is.na(diff.expressed.genes$name), ]
# diff.expressed.genes.map = queryMany(diff.expressed.genes$name, scopes="symbol", fields="entrezgene", species="human")
diff.expressed.genes.map = AnnotationDbi::select(EnsDb.Hsapiens.v86, keys = as.character(diff.expressed.genes$gene), keytype = "GENEID", columns = "ENTREZID")
# get the promoters of differentially expressed genes
diff.expressed.gene.promoters = promoters(txdb.genes[txdb.genes$gene_id %in% diff.expressed.genes.map$ENTREZID,])

# check to see if FXN (Entrez ID 2395) was in these promoters
diff.expressed.gene.promoters[which(diff.expressed.gene.promoters$gene_id == "2395")]

# get the percentage of overlap between the promoters (genes) of differentially expressed genes and
# the epigenetically marked (H3K4me3) promoters (genes) in the fetal brain cell line
# subsetByOverlaps() extracts the elements in the query (the first argument) that overlap at least one element in the subject (the second).
fetal.brain.overlap.H3K4me3 = subsetByOverlaps(diff.expressed.gene.promoters, fetal.brain.gr)
(fetal.brain.overlap.percentage.H3K4me3 = length(fetal.brain.overlap.H3K4me3) / length(diff.expressed.genes.map$ENTREZID) * 100)

# get the percentage of overlap between the promoters (genes) of differentially expressed genes and
# the epigenetically marked (H3K4me3) promoters (genes) in the adult brain cell line
# subsetByOverlaps() extracts the elements in the query (the first argument) that overlap at least one element in the subject (the second).
adult.brain.overlap.H3K4me3 = subsetByOverlaps(diff.expressed.gene.promoters, adult.brain.gr)
(adult.brain.overlap.percentage.H3K4me3 = length(adult.brain.overlap.H3K4me3) / length(diff.expressed.genes.map$ENTREZID) * 100)

# get the percentage of overlap between the promoters (genes) of differentially expressed genes and
# the epigenetically marked (H3K4me3) promoters (genes) in the liver cell line
# subsetByOverlaps() extracts the elements in the query (the first argument) that overlap at least one element in the subject (the second).
liver.line.overlap.H3K4me3 = subsetByOverlaps(diff.expressed.gene.promoters, liver.line.gr)
(liver.line.overlap.percentage.H3K4me3 = length(liver.line.overlap.H3K4me3) / length(diff.expressed.genes.map$ENTREZID) * 100)
      
