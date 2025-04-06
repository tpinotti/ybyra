# yplace: Y-chromosome phylogeny placement tool

yplace is a Snakemake workflow which calls Y-chromosome haplogroups by using a tree-based scoring of derived and ancestral SNP calls.



## Requirements

To use yplace you need:

- `python 3` (with `pandas`)
- `snakemake`
- `bcftools`
- A BAM file aligned to GRCh38

The Y-SNP tree used by yplace was built using yFull public data from June 2024 (v12.00.0), ensuring strict treeness for all informative SNPs.

## Getting Started

To run yplace, copy the `src/` and `tree/` folders into your working directory, and prepare `config.yml` and `units.tsv` files.

### `config.yml`

```
prefix:  # Project code name
ref: /path/to/hg38/genome.fa  # Path to reference genome
hg_snps: tree/jun24.snpinfo.yFull.bed
hg_info: tree/jun24.snpinfo.yFull.tsv
tree: tree/jun24_yFull.v12.tree
MQ: 30
units: units.tsv  # Path to units file
```

### `units.tsv`

```
sampleId    bam
```

- `sampleId`: the individual name  
- `bam`: path to the BAM file

### Running yplace

Once everything is set up, you can run the workflow, for example using 12 threads, like this:

```
snakemake -s chrY_v2.smk --configfile config.yml --cores 12
```


## Genotype Calling and Ancient DNA Damage

Genotypes are called using `bcftools`, requiring:

- At least 2 reads
- 70% majority to call a variant at any given locus

SNPs potentially affected by ancient DNA damage are flagged:

- 5' C>T
- 3' C>T (single-stranded library)
- 3' G>A (double-stranded library)

Library type can be specified (`ss`, `ds`, or `both`, with `both` as default). These flags do not affect placement by default, but a hard filter can be applied.



## Haplogroup Placement

For each node in the tree where a sample has derived SNPs, yplace calculates a tree score:

- +1 for every derived SNP from the root to the node
- –1 for every ancestral SNP along the same path

Instead of returning the node with the highest score, yplace selects the *optimal placement* — the highest-scoring node with no ancestral SNPs (i.e., 100% concordance).

If a downstream node has a higher score but includes ancestral calls, the sample is flagged as `unstable_downstream`.

### Score Ties

When multiple nodes have the same top tree score (often in low-coverage or low-resolution areas):

- The sample is flagged with `score_tie`
- yplace selects the node with the shortest path to the root
- If one of the tied nodes is ancestral to the others, it is also reported

---

## Main Output Files

### `aggregate.yplace`

This is the summary output table. Columns:

| Column | Description |
|--------|-------------|
| `individual` | Sample ID |
| `optplacement` | Selected placement node |
| `tree_score` | Total score |
| `flag` | Any flags set for this sample |
| `tree_path` | Tip-to-root path (root = `ybyra`) |

#### Flags

- `unstable_downstream`: Higher scoring node exists but with ancestral calls
- `score_tie`: Multiple nodes tied for top score
- `score_tie;shortest_path_to_root`: Selected node is the shortest path to root among ties
- `score_tie;shortest_path_to_root;most_recent_common_parent`: Also the shared ancestor of all tied nodes



### `unstabledownstream.yplace`

Lists nodes with a higher score than the optimal placement but containing ancestral SNPs.

| Column | Description |
|--------|-------------|
| `individual` | Sample ID |
| `id` | Node name |
| `derived` | Derived SNPs at node |
| `ancestral` | Ancestral SNPs at node |
| `tree_score` | Total score |
| `tree_path` | Tip-to-root path |



### `scoreties.yplace`

Lists all nodes with tied maximum scores.

| Column | Description |
|--------|-------------|
| `individual` | Sample ID |
| `id` | Node name |
| `derived` | Derived SNPs |
| `ancestral` | Ancestral SNPs |
| `tree_score` | Total score |
| `tree_path` | Tip-to-root path |



### `scoreties_summary.yplace`

Summarizes tie-breaking information per individual.

| Column | Description |
|--------|-------------|
| `individual` | Sample ID |
| `shortest_path_to_root` | Closest tied node to the root |
| `most_recent_common_parent` | Shared ancestor of tied nodes |

---

## Per Sample Output Files

Additional outputs per sample are generated across different folders:

### In `vcf/`

- `.vcf`: Informative Y-chromosome positions

### In `tables/`

- `.alleles.gz`: Compressed table with used VCF data

### In `calls/`

- `.calls`: All derived and ancestral SNP calls
- `.nopass`: Calls that failed filters

### In `yplace/`

- `.yplace`: Tree scores for all nodes in the phylogeny with a derived or ancestral call

---

## Acknowledgement

Thanks to Lucas Czech and Martin Sikora for bits and pieces of both code and ideas. All implementation problems, both of code and ideas, are mine. We also thank the yFull team for making their data freely available for the community.

## Citation

If you find this useful, you can cite yplace as:


https://www.biorxiv.org/content/10.1101/2024.03.13.584607v2




