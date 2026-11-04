#!/bin/bash
# Build transcript-to-gene table from the Ensembl cDNA FASTA headers
(echo "TXNAME,GENEID"
 zcat data/ref/Homo_sapiens.GRCh38.cdna.all.fa.gz | grep '^>' | sed 's/^>//' | \
 awk '{tx=$1; g="";
       for(i=1;i<=NF;i++){ if($i ~ /^gene:/){ split($i,a,":"); g=a[2] } }
       sub(/\.[0-9]+$/,"",tx); sub(/\.[0-9]+$/,"",g);
       print tx "," g}') > data/tx2gene.csv
