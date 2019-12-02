#!/bin/bash

# FETAL SAMPLES

echo "Aligning the reads of run SRR2071348"
hisat2 -p 8  --dta -x /data/hg38/genome -1 /data/sra/SRR2071348_1.fastq.gz -2 /data/sra/SRR2071348_2.fastq.gz \
		| samtools sort -@ 8 -n -o /data/align/SRR2071348.bam
echo "Finished the aignment of reads of run SRR2071348"

echo "Aligning the reads of run SRR2071349"
hisat2 -p 8  --dta -x /data/hg38/genome -1 /data/sra/SRR2071349_1.fastq.gz -2 /data/sra/SRR2071349_2.fastq.gz \
		| samtools sort -@ 8 -n -o /data/align/SRR2071349.bam
echo "Finished the aignment of reads of run SRR2071349"

echo "Aligning the reads of run SRR2071352"
hisat2 -p 8  --dta -x /data/hg38/genome -1 /data/sra/SRR2071352_1.fastq.gz -2 /data/sra/SRR2071352_2.fastq.gz \
		| samtools sort -@ 8 -n -o /data/align/SRR2071352.bam
echo "Finished the aignment of reads of run SRR2071352"

# ADULT SAMPLES

echo "Aligning the reads of run SRR2071346"
hisat2 -p 8  --dta -x /data/hg38/genome -1 /data/sra/SRR2071346_1.fastq.gz -2 /data/sra/SRR2071346_2.fastq.gz \
		| samtools sort -@ 8 -n -o /data/align/SRR2071346.bam
echo "Finished the aignment of reads of run SRR2071346"

echo "Aligning the reads of r SRR2071347"
hisat2 -p 8  --dta -x /data/hg38/genome -1 /data/sra/SRR2071347_1.fastq.gz -2 /data/sra/SRR2071347_2.fastq.gz \
		| samtools sort -@ 8 -n -o /data/align/SRR2071347.bam
echo "Finished the aignment of reads of run SRR2071347"

echo "Aligning the reads of run SRR2071350"
hisat2 -p 8  --dta -x /data/hg38/genome -1 /data/sra/SRR2071350_1.fastq.gz -2 /data/sra/SRR2071350_2.fastq.gz \
		| samtools sort -@ 8 -n -o /data/align/SRR2071350.bam
echo "Finished the aignment of reads of run SRR2071350"
