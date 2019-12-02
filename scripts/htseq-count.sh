#!/bin/bash

ANNOT=/data/gencode.v28.annotation.gtf

# FETAL SAMPLES

printf "\nn***** Assembling the transcripts for run SRR2071348 ...\n\n" 
htseq-count -f bam -s no -r name -m union --nonunique all /data/align/SRR2071348.bam $ANNOT > htseq/SRR2071348_counts.txt

printf "\nn***** Assembling the transcripts for run SRR2071349 ...\n\n" 
htseq-count -f bam -s no -r name -m union --nonunique all /data/align/SRR2071349.bam $ANNOT > htseq/SRR2071349_counts.txt

printf "\nn***** Assembling the transcripts for run SRR2071352 ...\n\n" 
htseq-count -f bam -s no -r name -m union --nonunique all /data/align/SRR2071352.bam $ANNOT > htseq/SRR2071352_counts.txt

# ADULT SAMPLES

printf "\nn***** Assembling the transcripts for run SRR2071346 ...\n\n" 
htseq-count -f bam -s no -r name -m union --nonunique all /data/align/SRR2071346.bam $ANNOT > htseq/SRR2071346_counts.txt

printf "\nn***** Assembling the transcripts for run SRR2071347  ...\n\n" 
htseq-count -f bam -s no -r name -m union --nonunique all /data/align/SRR2071347.bam $ANNOT > htseq/SRR2071347_counts.txt

printf "\nn***** Assembling the transcripts for run SRR2071350 ...\n\n" 
htseq-count -f bam -s no -r name -m union --nonunique all /data/align/SRR2071350.bam $ANNOT > htseq/SRR2071350_counts.txt


