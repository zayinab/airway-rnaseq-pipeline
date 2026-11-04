#!/bin/bash
mkdir -p results/salmon
for id in $(cat data/samples.txt); do
  echo "Quantifying $id ..."
  salmon quant -i data/ref/salmon_index -l A \
    -1 data/trimmed/${id}_1.fq.gz -2 data/trimmed/${id}_2.fq.gz \
    -p 4 --validateMappings -o results/salmon/${id}
done
