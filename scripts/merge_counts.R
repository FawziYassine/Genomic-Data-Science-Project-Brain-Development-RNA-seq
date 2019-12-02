# Set the working directory
setwd("~/brain/brain-zip/brain-zip/htseq")

# Load fetal samples
SRX683795 <- read.table("SRR2071348_counts.txt", header = FALSE)
SRX683796 <- read.table("SRR2071349_counts.txt", header = FALSE)
SRX683799 <- read.table("SRR2071352_counts.txt", header = FALSE)

# Load adult samples
SRX683793 <- read.table("SRR2071346_counts.txt", header = FALSE)
SRX683794 <- read.table("SRR2071347_counts.txt", header = FALSE)
SRX683797 <- read.table("SRR2071350_counts.txt", header = FALSE)

# check to see if all elements of the first column (transcripts) are the same across all 6 samples
all(SRX683797[,1] == SRX683794[,1])
all(SRX683793[,1] == SRX683794[,1])
all(SRX683799[,1] == SRX683793[,1])
all(SRX683797[,1] == SRX683794[,1])
all(SRX683795[,1] == SRX683797[,1])

# create a merged_counts table
merged_counts <- data.frame(row.names = SRX683795[,1], SRX683795 = SRX683795[,2], SRX683796 = SRX683796[,2]
                            , SRX683799 = SRX683799[,2], SRX683793 = SRX683793[,2]
                            , SRX683794 = SRX683794[,2], SRX683797 = SRX683797[,2])
write.table(merged_counts, "../../../merged_counts.tsv", quote = FALSE, sep = '\t')
