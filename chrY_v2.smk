import pandas as pd
from itertools import combinations
import glob
import os

## --------------------------------------------------------------------------------
## global parameters from config file
configfile: "config.yml"
PREFIX = config["prefix"] ## prefix for output files
REF = config["ref"] ## path to reference genomes in .fasta format
REF_BUILD = config.get("ref_build", "hg38") ## reference build: hg38, hg19, or hs37d5
MQ = config["MQ"] ## minimum mapping quality
HG_SNPS = config["hg_snps"] ## BED file with haplogroup SNP positions
HG_INFO = config["hg_info"] ## haplogroup info table
TREE = config["tree"] ## path to the tree file

# Liftover parameters (only needed if not hg38)
if REF_BUILD != "hg38":
    HG38_REF = config["hg38_ref"] ## path to hg38 reference
    CHAIN_FILE = config["chain_file"] ## path to hg19ToHg38 chain file (works for both hg19 and hs37d5)
    
## --------------------------------------------------------------------------------
## helpers
unit_df = pd.read_table(config["units"], comment="#").set_index(["sampleId"])
SAMPLES = unit_df.index.values.tolist()

# Determine chromosome name based on reference build
def get_chr_name():
    if REF_BUILD == "hs37d5":
        return "Y"
    else:  # hg19 or hg38
        return "chrY"

CHR_NAME = get_chr_name()

## --------------------------------------------------------------------------------
## functions
def get_bam(wildcards):
    return unit_df.loc[(wildcards.sample), "bam"]

def get_vcf_input(wildcards):
    """Return appropriate VCF input based on reference build"""
    if REF_BUILD == "hg38":
        return f"vcf/Y.{wildcards.sample}.{PREFIX}.vcf.gz"
    else:
        return f"vcf/Y.{wildcards.sample}.{PREFIX}.hg38.vcf.gz"

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
        f"samtools view -q{MQ} -F 1028 -bh {{input}} {CHR_NAME} > {{output}}"

rule get_vcf:
    input:
        get_bam
    output:
        "vcf/Y.{sample}." + PREFIX + ".vcf.gz"
    params:
        sm="{sample}"
    shell:
        f"""
        echo -e {{params.sm}} > vcf/{{params.sm}}.samples.txt
        echo -e '{{params.sm}}\t1' > vcf/{{params.sm}}.ploidy.txt
        bcftools mpileup -B -q {MQ} -f {REF} -r {CHR_NAME} --ignore-RG {{input}} | bcftools reheader -s vcf/{{params.sm}}.samples.txt | bcftools call -Am -Oz -S vcf/{{params.sm}}.ploidy.txt --ploidy GRCh38 > {{output}}
        """

rule index_vcf:
    input:
        "vcf/Y.{sample}." + PREFIX + ".vcf.gz"
    output:
        "vcf/Y.{sample}." + PREFIX + ".vcf.gz.tbi"
    shell:
        "bcftools index -t {input}"

# Liftover rule - only runs if reference is not hg38
rule liftover_vcf:
    input:
        vcf="vcf/Y.{sample}." + PREFIX + ".vcf.gz",
        tbi="vcf/Y.{sample}." + PREFIX + ".vcf.gz.tbi"
    output:
        vcf="vcf/Y.{sample}." + PREFIX + ".hg38.vcf.gz",
        reject="vcf/Y.{sample}." + PREFIX + ".reject.bcf"
    params:
        sm="{sample}"
    shell:
        f"""
        bcftools +liftover --no-version -Ou {{input.vcf}} -- \
          -s {REF} \
          -f {HG38_REF} \
          -c {CHAIN_FILE} \
          --reject {{output.reject}} \
          --reject-type b \
          --write-src \
          --drop-tags FORMAT/PL,FORMAT/AD | \
        bcftools sort -o {{output.vcf}} -Oz --write-index
        """

rule index_hg38_vcf:
    input:
        "vcf/Y.{sample}." + PREFIX + ".hg38.vcf.gz"
    output:
        "vcf/Y.{sample}." + PREFIX + ".hg38.vcf.gz.tbi"
    shell:
        "bcftools index -t {input}"

rule get_alleles:
    input:
        vcf=lambda wildcards: get_vcf_input(wildcards),
        tbi=lambda wildcards: get_vcf_input(wildcards) + ".tbi"
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

# Add conditional dependencies based on reference build
if REF_BUILD != "hg38":
    # If not hg38, we need liftover before get_alleles
    ruleorder: liftover_vcf > get_alleles