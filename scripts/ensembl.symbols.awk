


#!/bin/bash

awk -F "\t" '$3 == "gene" {print $9}' /data/Homo_sapiens.GRCh38.103.gtf | \
tr -d ";\"" | \
awk -F " " '{print $2, $6}'


