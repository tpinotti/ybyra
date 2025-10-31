# ybyra: Y-chromosome phylogeny placement tool

ybyra is a Snakemake workflow which calls Y-chromosome haplogroups from bam files by using a tree-based scoring of derived and ancestral SNP calls.

## Quick overview

With ybyra you can:

- Call Y-chromosome haplogroups from BAM files mapped to either hg37 or hg38
- Use either ISOGG, yFull or FamilyTreeDNA (FTDNA) Y-SNP trees
- Apply an optional ancient DNA (aDNA) damage filter
- Plot phylogenetic placements for all samples

## Requirements

To use ybyra you need:

- `python 3` (with `pandas`)
- `snakemake`
- `bcftools`
- `ete3` (only for plotting)
- A BAM file


## Getting Started

### 1. Get ybyra
To run ybyra, clone this repository or copy the `src/`, `trees/` and `configs/` folders into your working directory, as well as `ybyra.v03.smk`.

### 2. Create your  `units.tsv`
Prepare a tab-separated `units.tsv` file, listing your samples and bam paths:


### `units.tsv`

```
sampleId    bam
Kennewick	/projects/bam/rasmussen2015/Kennewick.bam	
```

- `sampleId`: individual name  
- `bam`: path to the BAM file


### 3. Choose and edit your configuration file
Then, we'll need to tell ybyra which Y-SNP tree we want to use, as well was which reference genome BAM files are mapped to. This is done through the yaml files in the `configs/` folder.


The `configs/` folder looks like this:

- `configs/hg37/hg37.isogg.yml`:  ISOGG Y-SNP tree in hg37
- `configs/hg38/hg38.isogg.yml`:  ISOGG Y-SNP tree in hg38
- `configs/hg37/hg37.yfull.yml`:  yFull Y-SNP tree in hg37
- `configs/hg38/hg38.yfull.yml`:  yFull Y-SNP tree in hg38
- `configs/hg37/hg37.ftdna.yml`: FamilyTree DNA Y-SNP tree in hg37
- `configs/hg38/hg38.ftdna.yml`: FamilyTree DNA Y-SNP tree in hg38

ybyra offer users three different Y-SNPs tree topologies (ISOGG, yFull and FamilyTreeDNA), which are based on both public and private datasets. As those datasets do not overlap, the tree topology is different. ybyra only uses SNPs occurring inside the 10Mb region defined in Poznik et al. 2013, and ensures strict treeness for all markers. Information for tree topology and SNPs (included or excluded after filtering) can be found in the `trees/` folder.


In this example, our bams are mapped to hg37 and we would like to use Family Tree DNA Y-SNP tree so we use `configs/hg37.ftdna.yml`. 


We then need to update the corresponding line in the config with the path to the human reference genome hg37 in your system. If your reference genome for some reason does not use Y (or chrY for hg38) to denote the Y-chromosome, a quick workaround is to modify the 'chrom' line accordingly.


### `configs/hg37.ftdna.yml`

```
prefix: ybyra
ref: path/to/reference/genome	## change here
chrom: Y
hg_snps: trees/ftdna/hg37/ftdnaY.oct25.hg37.bed
hg_info: trees/ftdna/hg37/ftdnaY.oct25.hg37.tsv
tree: trees/ftdna/ftdnaY.oct25.complete.tree
MQ: 30
units:  units.tsv
damage_filter: false
lib_type: both
```

### 4. (Optional) Enable Ancient DNA damage filter
Finally, ybyra natively has an ancient DNA damage filter. This again is done by editing your chosen config file. If you set `damage_filter` to `true`, ybyra will call haplogroups excluding all SNPs flagged as potentially deriving from ancient DNA damage.


As the damage profile is dependent on library type, users must report whether libraries were `ds` (double-stranded), `ss` (single-stranded) or `both` (both library types in the bam file or unknown; default).


If `damage_filter` is set to `false`, ybyra will still flag SNPs as damaged, but will not perform any filtering.


### 5. Run ybyra

Once everything is set up, you can run the workflow, for example, using 12 threads, like this:

```
snakemake -s ybyra.v03.smk --configfile configs/hg37/hg37.ftdna.yml --cores 12
```


## Genotype Calling and Ancient DNA Damage

Genotypes are called using `bcftools`, requiring 70% majority to call a variant at any given locus.

SNPs potentially affected by ancient DNA damage are flagged, following library type damage profile and read orientation.

- 5' C>T (forward – all libraries types)
- 3' C>T (reverse – only single-stranded libraries)
- 3' G>A (reverse – only double-stranded libraries)

If a SNP is still supported by a 70% majority after excluding the support from reads potentially affected by damaged, it is not flagged.


## Haplogroup Placement

For each node in the tree where a sample has a derived or ancestral SNPs, ybyra calculates a tree score:

- +1 for every derived SNP from the root to the node
- –1 for every ancestral SNP along the same path

Instead of just returning the node with the highest score, ybyra selects the *optimal placement* — the highest-scoring node which is supported by at least another derived hit 5 nodes upstream (what we call "5 step rule").

In the case of the optimal placement including ancestral calls, the sample is flagged as `unstable_downstream`.

### Score Ties

When multiple nodes have the same top tree score  (often in low-coverage or low-resolution areas) and all pass the 5 step rule:

- The sample is flagged as `score_tie`
- ybyra selects the Most Recent Common Ancestor (MRCA) of all tied nodes with a derived hit as optimal placement


## Plotting

ybyra generates two plots:

- A tree showing the optimal placement for each individual
- A second tree showing all tied-score placements

Example plots from ancient individuals from Antonio et al. 2019 (10.1126/science.aay6826) are in the `examples/` folder.

## Main Output Files

### `aggregate.yplace`

This is the summary output table. Columns:

| Column | Description |
|--------|-------------|
| `individual` | Sample ID |
| `optplacement` | Selected placement node |
| `tree_score` | Total score |
| `flag` | Any flags set for this sample |
| `tree_path` | Tip-to-root path (root = `ybyra`, 'tree' in Tupi) |

#### Flags

- `unstable_downstream`: Ancestral calls at optimal placement is different from 0
- `score_tie;most_recent_common_parent`: Multiple nodes tied for high score; optimal placement is the most recent common ancestor of the tied nodes with a derived hit
-  `tree_score_below_50`:  Low confidence placement (low-coverage)
- `step_rule`:  Higher scoring nodes failed the 5 step rule

If you prefer the `tree_path` direction to be root-to-tip instead, a little helper script `src/ynvert.py` takes an `aggregate.yplace` file and outputs it with inverted `tree_path`.   

### `aggregate.pdf`

This is the summary output plot, showing the full path and relationship between all individuals in `aggregate.yplace`.

Individuals with the `tree_score_below_50` flag are denoted with a ** symbol.

### `unstabledownstream.yplace`

Details on nodes of individuals where the optimal placement contains ancestral SNPs (flag: `unstable_downstream`).

| Column | Description |
|--------|-------------|
| `individual` | Sample ID |
| `id` | Node name |
| `derived` | Derived SNPs at node |
| `ancestral` | Ancestral SNPs at node |
| `tree_score` | Total score |
| `tree_path` | Tip-to-root path |


### `scoreties.yplace`

Lists all nodes with tied maximum scores (flag: `score_tie`).

| Column | Description |
|--------|-------------|
| `individual` | Sample ID |
| `id` | Node name |
| `derived` | Derived SNPs |
| `ancestral` | Ancestral SNPs |
| `tree_score` | Total score |
| `tree_path` | Tip-to-root path |


### `scoreties_summary.yplace`

Summarizes score ties information per individual.

| Column | Description |
|--------|-------------|
| `individual` | Sample ID |
| `shortest_path_to_root` | Closest tied node to the root |
| `most_recent_common_parent` | Shared ancestor of tied nodes |

### `scoreties.pdf`

This is the score ties summary output plot, showing the full path and relationship between all possible placement nodes of individuals in `scoreties.yplace`.


### `step5nopass.yplace`

List all nodes for all samples which failed the 5 step rule.

| Column | Description |
|--------|-------------|
| `individual` | Sample ID |
| `id` | Node name |
| `tree_score` | Total score |
| `tree_path` | Tip-to-root path |

### `fail.yplace`

List all individuals with `tree_score` below 10.

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


## Other useful files

The `trees/` directory includes complete tree structures for the three topologies, along with SNP lists and their node placements.

SNPs excluded from analysis (non-unique or outside the 10 Mb region) are listed in `trees/*/hg38/discarded/`.

Liftover from hg38 to hg37 was performed using CrossMap (https://github.com/liguowang/CrossMap). Code is available at`liftover/liftoverY.py`. SNPs failing liftover are listed under  `trees/*/hg37/liftoverfail/`.


## Acknowledgement

Thanks to Lucas Czech and Martin Sikora for code and ideas, and to Hugh McColl, Teemu and Armando for helpful suggestions and comments.

 We also thank ISOGG (https://isogg.org),  yFull (https://www.yfull.com) and FamilyTree DNA (https://discover.familytreedna.com) for making their tree publicly available for the community.

This is a stable draft version of ybyra; ideas, suggestions and comments are very welcome. You can get in touch at thomaz.pinotti(at)sund.ku.dk

## Version

This is ybyra v0.3 (30/10/25).

## Citation

If you find this useful, for now you can cite ybyra as:

https://www.biorxiv.org/content/10.1101/2024.03.13.584607v2




