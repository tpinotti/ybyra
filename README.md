# ybyra: Y-chromosome haplogroup calling using a tree-based scoring method

[![License](https://img.shields.io/github/license/tpinotti/ybyra)](https://github.com/tpinotti/ybyra/blob/main/LICENSE.md)
[![CI](https://github.com/tpinotti/ybyra/actions/workflows/ci.yaml/badge.svg)](https://github.com/tpinotti/ybyra/actions/workflows/ci.yaml)
[![Release](https://img.shields.io/github/v/release/tpinotti/ybyra)](https://github.com/tpinotti/ybyra/releases)
[![DOI](https://img.shields.io/badge/DOI-10.1101%2F2025.11.20.689455-blue)](https://doi.org/10.1101/2025.11.20.689455)

ybyra is a Snakemake workflow which calls Y-chromosome haplogroups from bam files by using a tree-based scoring of derived and ancestral SNP calls.


## Quick overview

With ybyra you can:

- Call Y-chromosome haplogroups from BAM files mapped to either hg37 or hg38
- Use either ISOGG, YFull or FamilyTreeDNA (FTDNA) Y-SNP trees
- Apply an optional ancient DNA (aDNA) damage filter
- Plot phylogenetic placements for all samples


## Requirements

To use ybyra you need:

- `python <= 3.12` (with `pandas`)
- `snakemake`
- `bcftools`
- `ete3 <= 3.1.3` (only for plotting)

A conda env file with these dependencies is provided in `workflow/envs/ybyra.yaml`. Use for instance `micromamba create -f workflow/envs/ybyra.yaml -n ybyra` (or your preferred conda-like tool) to create an environment, and `micromamba activate ybyra` to activate it prior to running the steps below.


## Getting Started

For a complete, ready-to-run example with ancient DNA samples, see the [`example/`](https://github.com/tpinotti/ybyra/tree/main/example) folder.

### 1. Get ybyra

To run ybyra, clone this repository or download a stable [release version](https://github.com/tpinotti/ybyra/releases) or the latest update from the green "Code" button above via ["Download ZIP"](https://github.com/tpinotti/ybyra/archive/refs/heads/main.zip).


### 2. Create a working directory

First, we create a *working directory* where we want our analysis to run, i.e., where output files will be generated. It can be called anything, such as `my_analysis`.


### 3. Create your `units.tsv`

Prepare a tab-separated `units.tsv` file in the working directory, listing your samples and bam paths:

### `units.tsv`

```
sampleId    bam
Kennewick	/projects/bam/rasmussen2015/Kennewick.bam
```

- `sampleId`: individual name
- `bam`: path to the BAM file

Note that bam files need to be indexed.

### 4. Create your configuration file

Then, we'll need to tell ybyra which Y-SNP tree we want to use, as well was which reference genome the BAM files are mapped to. This is done through the config file in `config/config.yaml`. Copy this file to the working directory, keeping the name `config.yaml`. Then, edit it as needed - explanations of all settings are give in the file itself.

Alternatively, the `config.yaml` in the working directory can contain only the settings that differ from the defaults in `config/config.yaml`; all other settings are taken from there. Unknown settings (e.g., typos or settings from outdated config files) are reported as an error.

In particular, the path to the reference genome, as well as its build and the reference SNP tree need to be changed as needed. ybyra offer users three different Y-SNPs tree topologies (ISOGG, YFull and FamilyTreeDNA), which are based on both public and private datasets. As those datasets do not overlap, the tree topology is different. ybyra only uses SNPs occurring inside the 10Mb region defined in [Poznik et al. 2013](https://doi.org/10.1126/science.1237619), and ensures strict treeness for all markers. Information for tree topology and SNPs (included or excluded after filtering) can be found in the `trees/` folder.

The three available tree topologies are distributed under different licenses, all of which require proper attribution. Users are responsible for ensuring that the appropriate source is credited when a tree is used. Details on licensing and attribution can be found in the [attribution](https://github.com/tpinotti/ybyra/tree/main#attribution) section of this documentation and in the [`trees`](https://github.com/tpinotti/ybyra/tree/main/trees) subdirectory.

In this example, our bams are mapped to hg37 and we would like to use FamilyTreeDNA Y-SNP tree so we use `build: "hg37"` and `tree: "ftdna"` in the config file, and need to update the line in the config with the path to the human reference genome hg37 in your system.


### 5. (Optional) Enable Ancient DNA damage filter

Finally, ybyra natively has an ancient DNA damage filter. This again is done by editing your config file. If you set `damage_filter` to `true`, ybyra will call haplogroups after excluding SNPs flagged as potentially deriving from ancient DNA damage.

As the damage profile is dependent on library type, users must report whether libraries were `ds` (double-stranded), `ss` (single-stranded) or `both` (both library types in the bam file or unknown; default).

Furthermore, when assessing damage, ybyra supports three different damage models. Their difference lies on how ancestral calls are treated. The damage models can be chosen in the config (`damage_model`):

- `naive` – considers all C>T (and G>A depending on library type) as damaged
- `uni` – unidirecional damage model. considers C>T (and G>A depending on library type) as damaged if derived; ancestral calls are never considered damaged
- `bi` – bidirecional damage model. considers C>T (and G>A depending on library type) as damaged if derived; instead consider T>C (and A>G depending on library type) as damaged if ancestral

If `damage_filter` is set to `false`, ybyra will still flag SNPs as damaged, but will not perform any filtering.


### 6. Run ybyra

At this stage, the working directory should contain two files, `config.yaml` and `units.tsv`.
Once everything is set up, you can run the workflow from the ybyra directory, for example, using 12 threads, like this:

```
cd /path/to/ybyra
snakemake --cores 12 --directory /path/to/my_analysis
```

That is, we run ybyra while in the directory with the code, and then use the `--directory` option of snakemake to point to the working directory (here assuming `my_analysis` as above, but replace as needed). This is where result files will then be produced, as explained in the following.


## Genotype Calling and Ancient DNA Damage

Genotypes are called using `bcftools`, requiring at least `min_majority_fraction` (default: 70%) of the reads to support either the reference or the alternative allele at any given locus.

SNPs potentially affected by ancient DNA damage are flagged, following library type and damage model:

- `damage_filter: false` : uses all SNPs in scoring

- `damage_filter: true`; +  `lib_type`: `"both"` or `"dslib"` + `damage_model: "naive"` : excludes from scoring all C>T and G>A SNPs, regardless if derived or ancestral

- `damage_filter: true` + `lib_type: "sslib"` + `damage_model: "naive"` : excludes from scoring all C>T SNPs, regardless if derived or ancestral

- `damage_filter: true` + `lib_type`: `"both"` or `"dslib"` + `damage_model: "uni"` : excludes from scoring all C>T and G>A derived SNPs

- `damage_filter: true` + `lib_type: "sslib"` + `damage_model: "uni"` : excludes from scoring all C>T derived SNPs

- `damage_filter: true` + `lib_type`: `"both"` or `"dslib"` + `damage_model: "bi"` : excludes from scoring all C>T and G>A derived SNPs and all T>C and A>G ancestral SNPs

- `damage_filter: true` + `lib_type: "sslib"` + `damage_model: "bi"` : excludes from scoring all C>T derived SNPs and all T>C ancestral SNPs


## Haplogroup Placement

For each node in the tree where a sample has a derived or ancestral SNPs, ybyra calculates a tree score:

- +1 for every derived SNP from the root to the node
- –1 for every ancestral SNP along the same path

Instead of just returning the node with the highest score, ybyra selects the *optimal placement* — the highest-scoring node which is supported by at least one other derived hit within `step_size` nodes upstream (default: 5), what we call the "step rule". Setting `step_size: 0` disables the rule.

In the case of the optimal placement including ancestral calls, the sample is flagged as `unstable_downstream`.


### Score Ties

When multiple nodes have the same top tree score  (often in low-coverage or low-resolution areas) and all pass the step rule:

- The sample is flagged as `score_tie`
- ybyra selects the Most Recent Common Ancestor (MRCA) of all tied nodes with a derived hit as optimal placement

### Exploring the Step Size

The step rule trades resolution for robustness: smaller step sizes avoid placements on wrong lineages in low-coverage samples, but place more samples further upstream. To see how your samples are affected, set `step_size_exploration: true`. This places all samples at step sizes 1 to `step_size_exploration_max` (default: 10) and without the step rule, and compares them to the placement at your `step_size` (or without the step rule, see `step_size_exploration_baseline`). Output in `step_size/`:

- `step_size.tsv`: Placement of each sample at each step size, with its class relative to the baseline (identical, upstream, downstream, other lineage, or failed) and the sample depth (average reads per panel SNP)
- `overview.pdf`: Fraction of samples per class for each step size, and how many placements change, by sample depth
- `samples.pdf`: Placements of all samples that change with the step size


## Plotting

ybyra generates two plots:

- A tree showing the optimal placement for each individual
- A second tree showing all tied-score placements

Example plots from ancient individuals from [Antonio et al. 2019](https://doi.org/10.1126/science.aay6826) are in the `gallery/` folder.

For publication-ready SVG trees of each sample and a summary of all samples, see our plotting tool [ybyrapûera](puera/).


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
- `low_tree_score`: Low confidence placement (low-coverage), with a tree score below `low_tree_score` (default: 50)
- `step_rule`:  Higher scoring nodes failed the step rule

If you prefer the `tree_path` direction to be root-to-tip instead, a little helper script `workflow/scripts/ynvert.py` takes an `aggregate.yplace` file and outputs it with inverted `tree_path`.

### `aggregate.pdf`

This is the summary output plot, showing the full path and relationship between all individuals in `aggregate.yplace`.

Individuals with the `low_tree_score` flag are denoted with a ** symbol.

### `unstable_downstream.yplace`

Details on nodes of individuals where the optimal placement contains ancestral SNPs (flag: `unstable_downstream`).

| Column | Description |
|--------|-------------|
| `individual` | Sample ID |
| `id` | Node name |
| `derived` | Derived SNPs at node |
| `ancestral` | Ancestral SNPs at node |
| `tree_score` | Total score |
| `tree_path` | Tip-to-root path |


### `score_ties.yplace`

Lists all nodes with tied maximum scores (flag: `score_tie`).

| Column | Description |
|--------|-------------|
| `individual` | Sample ID |
| `id` | Node name |
| `derived` | Derived SNPs |
| `ancestral` | Ancestral SNPs |
| `tree_score` | Total score |
| `tree_path` | Tip-to-root path |


### `score_ties_summary.yplace`

Summarizes score ties information per individual.

| Column | Description |
|--------|-------------|
| `individual` | Sample ID |
| `shortest_path_to_root` | Closest tied node to the root |
| `most_recent_common_parent` | Shared ancestor of tied nodes |

### `score_ties.pdf`

This is the score ties summary output plot, showing the full path and relationship between all possible placement nodes of individuals in `score_ties.yplace`.


### `step_rule_nopass.yplace`

List all nodes for all samples which failed the step rule. Empty (header only) when the rule is disabled.

| Column | Description |
|--------|-------------|
| `individual` | Sample ID |
| `id` | Node name |
| `tree_score` | Total score |
| `tree_path` | Tip-to-root path |

### `fail.yplace`

List all individuals with `tree_score` below `min_tree_score` (default: 10), with one of these flags:

- `below_min_tree_score`: No node reaches the minimum tree score
- `below_min_tree_score_after_step_rule`: All nodes above the minimum tree score failed the step rule

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

The `trees/` directory in the main code directory includes complete tree structures for the three topologies, along with SNP lists and their node placements.

SNPs excluded from analysis (non-unique or outside the 10 Mb region) are listed in `trees/*/hg38/discarded/`.

Liftover from hg38 to hg37 was performed using CrossMap (https://github.com/liguowang/CrossMap). Code is available at`trees/liftover/liftoverY.py`. SNPs failing liftover are listed under  `trees/*/hg37/liftoverfail/`.


## Acknowledgement

Thanks to J. Víctor Moreno Mayar, Teemu ([@teepean](https://github.com/teepean)) and Armando for helpful suggestions and comments.

We also thank ISOGG ([isogg.org](https://isogg.org)),  YFull ([yfull.com](https://www.yfull.com)) and FamilyTreeDNA ([discover.familytreedna.com](https://discover.familytreedna.com)) for making their trees available for the community.

Ideas, suggestions and comments are very welcome. You can get in touch at thomaz.pinotti(at)sund.ku.dk, or open an [Issue](https://github.com/tpinotti/ybyra/issues) here on GitHub.

## Attribution

ybyra is published under the [MIT License](https://github.com/tpinotti/ybyra/blob/dev/LICENSE.md).

ybyrapûera, in the [`puera`](puera/) subdirectory, is published under the [GPLv3](puera/LICENSE.txt).

International Society of Genetic Genealogy (ISOGG) tree is distributed under the [CC BY-NC-SA 3.0](https://creativecommons.org/licenses/by-nc-sa/3.0/deed.en) license.

YFull tree is distributed under the [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/deed.en) license.

FamilyTreeDNA Y-DNA Haplotree is distributed under the [CC BY-NC-ND 4.0](https://creativecommons.org/licenses/by-nc-nd/4.0/deed.en) license.

Details on the licenses and tree versions can be found in the [`trees`](https://github.com/tpinotti/ybyra/tree/main/trees) subdirectory README.


## Citation

If you find ybyra useful, for now you can cite its preprint:

> **ybyra: Y-chromosome haplogroup calling using a tree-based scoring method**.<br />
>  Thomaz Pinotti, Hugh McColl, Martin Sikora, Lucas Czech.<br />
> *bioRxiv*, DOI:[10.1101/2025.11.20.689455](https://doi.org/10.1101/2025.11.20.689455)
