SAMPLES = [line.strip() for line in open("data/samples.txt") if line.strip()]

rule all:
    input:
        "results/deseq2_results.csv",
        "results/figures/pca.png",
        "results/qc/raw/multiqc_raw.html"

rule fastqc_raw:
    input:
        expand("data/raw/{s}_{r}.fastq.gz", s=SAMPLES, r=[1, 2])
    output:
        "results/qc/raw/multiqc_raw.html"
    threads: 4
    shell:
        "bash scripts/02_qc.sh"

rule trim:
    input:
        r1="data/raw/{sample}_1.fastq.gz",
        r2="data/raw/{sample}_2.fastq.gz"
    output:
        r1="data/trimmed/{sample}_1.fq.gz",
        r2="data/trimmed/{sample}_2.fq.gz",
        html="results/qc/fastp/{sample}_fastp.html",
        json="results/qc/fastp/{sample}_fastp.json"
    threads: 4
    shell:
        "fastp -i {input.r1} -I {input.r2} -o {output.r1} -O {output.r2} "
        "-h {output.html} -j {output.json} -w {threads}"

rule salmon_quant:
    input:
        r1="data/trimmed/{sample}_1.fq.gz",
        r2="data/trimmed/{sample}_2.fq.gz",
        index="data/ref/salmon_index"
    output:
        "results/salmon/{sample}/quant.sf"
    threads: 4
    shell:
        "salmon quant -i {input.index} -l A -1 {input.r1} -2 {input.r2} "
        "-p {threads} --validateMappings -o results/salmon/{wildcards.sample}"

rule deseq2:
    input:
        quant=expand("results/salmon/{s}/quant.sf", s=SAMPLES),
        meta="metadata/metadata.csv",
        tx2gene="data/tx2gene.csv"
    output:
        "results/deseq2_results.csv",
        "results/figures/pca.png",
        "results/figures/volcano.png",
        "results/figures/heatmap_top30.png"
    shell:
        "Rscript scripts/07_deseq2.R"
