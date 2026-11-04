#!/bin/bash
# Trim adapters and low-quality bases with fastp
mkdir -p data/trimmed results/qc/fastp
for id in $(cat data/samples.txt); do
  echo "Trimming $id ..."
  fastp -i data/raw/${id}_1.fastq.gz -I data/raw/${id}_2.fastq.gz \
        -o data/trimmed/${id}_1.fq.gz -O data/trimmed/${id}_2.fq.gz \
        -h results/qc/fastp/${id}_fastp.html \
        -j results/qc/fastp/${id}_fastp.json -w 4
done
multiqc results/qc/fastp -o results/qc/fastp -n multiqc_trimmed
