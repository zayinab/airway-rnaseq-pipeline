#!/bin/bash
# Quality control of raw reads
mkdir -p results/qc/raw
fastqc data/raw/*.fastq.gz -o results/qc/raw -t 4
multiqc results/qc/raw -o results/qc/raw -n multiqc_raw
