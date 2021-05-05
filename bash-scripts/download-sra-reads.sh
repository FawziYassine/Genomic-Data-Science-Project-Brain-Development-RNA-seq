#!/bin/bash

# FETAL SAMPLES

fastq-dump -v --gzip --split-files -O /data/sra  SRR2071348
printf “\n****** SRR2071348 is downloaded  ...\n\n”

fastq-dump -v --gzip --split-files -O /data/sra  SRR2071349
printf “\n****** SRR2071349 is downloaded  ...\n\n”

fastq-dump -v --gzip --split-files -O /data/sra  SRR2071352
printf “\n****** SRR2071352 is downloaded  ...\n\n”

# ADULT SAMPLES

fastq-dump -v --gzip --split-files -O /data/sra  SRR2071346
printf “\n****** SRR2071346 is downloaded  ...\n\n”

fastq-dump -v --gzip --split-files -O /data/sra  SRR2071347
printf “\n****** SRR2071347 is downloaded  ...\n\n”

fastq-dump -v --gzip --split-files -O /data/sra SRR2071350
printf “\n****** SRR2071350 is downloaded  ...\n\n”

