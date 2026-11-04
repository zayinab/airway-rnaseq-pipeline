#!/bin/bash
salmon index -t data/ref/Homo_sapiens.GRCh38.cdna.all.fa.gz \
             -i data/ref/salmon_index -k 31 -p 4
