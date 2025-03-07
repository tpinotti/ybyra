import pandas as pd
from itertools import combinations
import glob

## --------------------------------------------------------------------------------
## global parameters from config file
configfile: "config.yml"
PREFIX = config["prefix"] ## prefix for output files
REF = config["ref"] ## path to reference genomes in .fasta format
MQ = config["MQ"] ## minimum mapping quality
HG_SNPS = config["hg_snps"] ## BED file with haplogroup SNP positions
HG_INFO = config["hg_info"] ## haplogroup info table
TREE = config["tree"] ## path to the tree file

## --------------------------------------------------------------------------------
## helpers

unit_df = pd.read_table(config["units"], comment="#").set_index(["sampleId"])
SAMPLES = unit_df.index.values.tolist()

## --------------------------------------------------------------------------------
## functions

def get_bam(wildcards):
    return unit_df.loc[(wildcards.sample), "bam"]

## --------------------------------------------------------------------------------
## output file sets

callfiles = expand("calls/{sample}.calls", sample=SAMPLES)
yplacefiles = expand("yplace/{sample}.yplace", sample=SAMPLES)
aggfile = "aggregate.yplace"

## --------------------------------------------------------------------------------
## targets

rule all:
    input:
        callfiles,
        yplacefiles,
        aggfile

## --------------------------------------------------------------------------------
## rules
rule get_y_bam:
    input:
        get_bam
    output:
        "bam/Y.{sample}." + PREFIX + ".bam"
    shell:
        "samtools view -q{MQ} -F 1028 -bh {input} chrY > {output}"

rule get_vcf:
    input:
        get_bam
    output:
        "vcf/Y.{sample}." + PREFIX + ".vcf.gz"
    params:
        sm="{sample}"
    shell:
        """
        echo -e {params.sm} > vcf/{params.sm}.samples.txt
        echo -e '{params.sm}\t1' > vcf/{params.sm}.ploidy.txt
        bcftools mpileup -B -q {MQ} -f {REF} -r chrY --ignore-RG {input} | bcftools reheader -s vcf/{params.sm}.samples.txt | bcftools call -Am -Oz -S vcf/{params.sm}.ploidy.txt --ploidy GRCh38 > {output}
        """

rule index_vcf:
    input:
        "vcf/Y.{sample}." + PREFIX + ".vcf.gz"
    output:
        "vcf/Y.{sample}." + PREFIX + ".vcf.gz.tbi"
    shell:
        """
        bcftools index -t {input}
        """

rule get_alleles:
    input:
        vcf="vcf/Y.{sample}." + PREFIX + ".vcf.gz",
        tbi="vcf/Y.{sample}." + PREFIX + ".vcf.gz.tbi",
    output:
        "tables/Y.{sample}." + PREFIX + ".alleles.gz"
    shell:
        """
        bcftools query -T {HG_SNPS} -e 'TYPE=="indel"' -f '%POS\t%REF\t%ALT\t%DP\t%QUAL\t%DP4\t[%TGT]\n' {input.vcf} | gzip > {output}
        """

rule Ytree:
     input:
        alleles="tables/Y.{sample}." + PREFIX + ".alleles.gz"
     output:
       calls="calls/{sample}.calls",
       nopass="calls/{sample}.nopass",
     params:
       sm="{sample}",
     shell:
       "python src/ytree.py --lib both --alleles {input} --out calls/{params.sm} --snpinfo {HG_INFO}"

rule Yplace:
    input:
        calls="calls/{sample}.calls"
    output:
        "yplace/{sample}.yplace"
    shell:
        "python src/yplace.py {TREE} {input.calls} {output}"

rule aggregate_yplace:
    input:
        yplacefiles
    output:
        aggfile
    shell:
        "python src/ysummary.py {input} > {output}"
