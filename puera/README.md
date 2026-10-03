ybyrapûera
==========

ybyrapûera (`puera`) plots the results of a ybyra run as SVG trees, using [genesis](https://github.com/lczech/genesis).

## Setup

You need CMake and a C++17 compiler. Then, in this directory, run `make`, which builds `bin/puera`.
The genesis submodule is used if checked out, and downloaded otherwise.

## Usage

`puera sample` plots one tree per sample, with branches colored by the tree score of the sample, and its placement marked by a ring (dashed for score ties). For instance, after running the [`example`](../example), from the main directory:

```
puera/bin/puera sample \
    --tree-table-file trees/ftdna/ftdnaY.oct25.complete.tree \
    --snps-table-file trees/ftdna/hg37/ftdnaY.oct25.hg37.tsv \
    --clade-labels-file trees/ftdna/ftdnaY.oct25.labels.tsv \
    --ybyra-dir example \
    --exclude-damage \
    --out-dir plots
```

`puera summary` plots where all samples are placed, with branches colored by the number of samples placed in their clade. With `--groups-file`, it also plots each group of samples:

```
puera/bin/puera summary \
    --tree-table-file trees/ftdna/ftdnaY.oct25.complete.tree \
    --snps-table-file trees/ftdna/hg37/ftdnaY.oct25.hg37.tsv \
    --clade-labels-file trees/ftdna/ftdnaY.oct25.labels.tsv \
    --ybyra-dir example \
    --out-dir plots
```

Use the tree of the ybyra run. The SNP table is optional, and used for branch lengths. `--exclude-damage` needs to match the `damage_filter` setting of the run. See `--help` for all options.

The `trees/*/*.labels.tsv` files label the major haplogroups. They are made with `trees/make_clade_labels.py`, and can be edited by hand.

## License

ybyrapûera is published under the [GPLv3](LICENSE.txt), as required by genesis. Note that this differs from the [MIT License](../LICENSE.md) of ybyra itself.
