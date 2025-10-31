import pandas as pd
from itertools import combinations
import glob

## --------------------------------------------------------------------------------
## global parameters from config file
PREFIX = config["prefix"] ## prefix for output files
REF = config["ref"] ## path to reference genomes in .fasta format
MQ = config["MQ"] ## minimum mapping quality
HG_SNPS = config["hg_snps"] ## BED file with haplogroup SNP positions
HG_INFO = config["hg_info"] ## haplogroup info table
TREE = config["tree"] ## path to the tree file
CHROM = config["chrom"]

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
aggplot = "aggregate.pdf"
scoretiesplot = "scoreties.pdf"

## --------------------------------------------------------------------------------
## targets

rule all:
    input:
        callfiles,
        yplacefiles,
        aggfile,
        aggplot,
        scoretiesplot

## --------------------------------------------------------------------------------
## rules
rule get_y_bam:
    input:
        get_bam
    output:
        "bam/Y.{sample}." + PREFIX + ".bam"
    shell:
        "samtools view -q{MQ} -F 4 -bh {input} {CHROM} > {output}"

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
        bcftools mpileup -B -q {MQ} -f {REF} -r {CHROM} --ignore-RG {input} | bcftools reheader -s vcf/{params.sm}.samples.txt | bcftools call -Am -Oz --ploidy 1 > {output}
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
       lib = lambda wildcards: config.get("lib_type", "both")
     shell:
       "python src/ytree.py --lib {params.lib} --alleles {input} --out calls/{params.sm} --snpinfo {HG_INFO}"

rule Yplace:
    input:
        calls="calls/{sample}.calls"
    output:
        "yplace/{sample}.yplace"
    params:
        damage = lambda wildcards: "--no-damage" if config.get("damage_filter", False) else ""
    shell:
        "python src/yplace.py {params.damage} {TREE} {input.calls} {output}"

rule aggregate_yplace:
    input:
        yplacefiles
    output:
        aggfile
    shell:
        "python src/ysummary.py {input} > {output}"

rule plot_aggregate:
    input:
        aggfile
    output:
        "aggregate.pdf"
    shell:
        "python src/yplot.py {input} {output}"

rule plot_scoreties:
    input:
        agg = aggfile
    output:
        "scoreties.pdf"
    shell:
        "python src/yplot.py scoreties.yplace {output}"


