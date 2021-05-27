configfile: "config.yaml"

import os
import glob

SAMPLES = "/Users/angelol/Documents/PhD/BE-VCF/Christos_2021/data/*"
SAMPLES = sorted([os.path.splitext(val)[0] for val in (glob.glob(SAMPLES))]) #Remove .gz from filename path
SAMPLES = [os.path.splitext(val)[0] for val in SAMPLES]
SAMPLES = [os.path.basename(val) for val in SAMPLES]

for i, s in enumerate(SAMPLES):
    SAMPLES[i] = SAMPLES[i][:-12] # Remove _1 or _2 suffix for paired-end reads

SAMPLES = list(set(SAMPLES))

rule all:
    input:
        "calls/mpileup.tsv"

rule qual_filter:
    input:
        R1 = config["paths"]["raw_reads"] + "{sample}_L001_R1_001.fastq.gz",
        R2 = config["paths"]["raw_reads"] + "{sample}_L001_R2_001.fastq.gz"
    output:
        R1 = "filtered_reads/" + "{sample}_L001_R1_001.fastq.gz",
        R2 = "filtered_reads/" + "{sample}_L001_R2_001.fastq.gz"
    log:
        "logs/fastp/{sample}.html"
    shell:
        "fastp --thread 8 --correction -q 30 "
        "--in1 {input.R1} --in2 {input.R2} "
        "--out1 {output.R1} --out2 {output.R2} "
        "-h {log}"
        
rule bwa_map:
    input:
        template = config["paths"]["template"] + "template.fa",
        R1 = "filtered_reads/" + "{sample}_L001_R1_001.fastq.gz",
        R2 = "filtered_reads/" + "{sample}_L001_R2_001.fastq.gz"
    output:
        "mapped_reads/{sample}.bam"
    shell:
        "bwa mem -t 4 {input.template} {input.R1} {input.R2} | "
        "samtools view -Sb - > {output}"

rule samtools_sort:
    input:
        "mapped_reads/{sample}.bam"
    output:
        "sorted_reads/{sample}.bam"
    shell:
        "samtools sort -T sorted_reads/{wildcards.sample} "
        "-O bam {input} > {output}"

rule samtools_index:
    input:
        "sorted_reads/{sample}.bam"
    output:
        "sorted_reads/{sample}.bam.bai"
    shell:
        "samtools index {input}"

rule samtools_mpileup:
    input:
        fa = config["paths"]["template"] + "template.fa",
        bam = expand("sorted_reads/{sample}.bam", sample=SAMPLES),
        bai = expand("sorted_reads/{sample}.bam.bai", sample=SAMPLES)
    output:
        "calls/mpileup.tsv"
    shell:
        "bcftools mpileup --max-depth 500000 -O v -a FORMAT/AD -f {input.fa} {input.bam} > {output}" # Setting max-depth to 0 allows INF reads at each genomic position.

#"freebayes --min-alternate-count 1 --min-alternate-fraction 0 -f {input.fa} --ploidy 1 {input.bam} > calls/all.vcf"
#"freebayes -f {input.fa} --ploidy 1 {input.bam} > all.vcf"
#"bcftools mpileup --max-depth 200000 -Ou -f {input.fa} {input.bam} | bcftools call --ploidy 1 -p 1 -cv - > {output}"
