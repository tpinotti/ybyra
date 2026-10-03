#!/usr/bin/env python3
"""
Make a clade labels file for `puera --clade-labels-file`, labelling the major Y haplogroups.

For each haplogroup letter B-T, the label goes to the node named exactly as the letter (ISOGG,
yFull), or otherwise to the topmost node named with the letter as prefix, e.g. `R-M207` (FTDNA).
The result is written to stdout, and meant to be checked and edited by hand if needed.

Usage: make_clade_labels.py <tree> > <labels.tsv>
"""

import csv
import string
import sys

LETTERS = string.ascii_uppercase[1:20]  # B-T; A is paraphyletic, and hence not labelled.
ROOT = "ybyra"


def main(tree_file):
    with open(tree_file) as f:
        parent = {row["id"]: row["parent"] for row in csv.DictReader(f)}

    # Depth of each node below the root; None for nodes that are not connected to it.
    depth = {}
    def set_depth(node):
        path = []
        while node not in depth and node != ROOT:
            if node not in parent or node in path:
                break
            path.append(node)
            node = parent[node]
        base = 0 if node == ROOT else depth.get(node)
        for i, n in enumerate(reversed(path)):
            depth[n] = None if base is None else base + i + 1
    for node in parent:
        set_depth(node)
    unconnected = sum(1 for d in depth.values() if d is None)
    if unconnected:
        print(f"Warning: {unconnected} nodes are not connected to the root, and ignored",
              file=sys.stderr)

    print("node\tlabel")
    for letter in LETTERS:
        if depth.get(letter) is not None:
            candidates = [letter]
        else:
            matches = [n for n in parent if n.startswith(letter + "-") and depth[n] is not None]
            top = min((depth[n] for n in matches), default=None)
            candidates = [n for n in matches if depth[n] == top]
        if len(candidates) != 1:
            print(f"{letter}: {len(candidates)} candidate nodes {candidates}, skipped",
                  file=sys.stderr)
            continue
        print(f"{candidates[0]}\t{letter}")
        print(f"{letter}: {candidates[0]} (parent {parent[candidates[0]]})", file=sys.stderr)


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit(__doc__.strip())
    main(sys.argv[1])
