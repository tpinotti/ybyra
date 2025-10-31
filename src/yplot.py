import sys
from ete3 import Tree, TreeStyle, NodeStyle, TextFace

input_file = sys.argv[1]
output_file = sys.argv[2]
tree_col = "tree_path"
sample_col = "individual"
flag_col = "flag"

paths = []
samples = []

with open(input_file) as f:
    header = f.readline().strip().split("\t")
    tree_idx = header.index(tree_col)
    sample_idx = header.index(sample_col)
    flag_idx = header.index(flag_col) if flag_col in header else None

    for line in f:
        fields = line.strip().split("\t")
        if len(fields) <= max(tree_idx, sample_idx):
            continue

        # invert tree path
        inverted_path = list(reversed(fields[tree_idx].split("<")))
        paths.append(inverted_path)

        # Append "**" if flag column exists and contains "tree_score_below_50"
        sample_name = fields[sample_idx]
        if flag_idx is not None and "tree_score_below_50" in fields[flag_idx]:
            sample_name += "**"
        samples.append(sample_name)


def build_tree(paths):
    root = {}
    for path in paths:
        node = root
        for p in path:
            if p not in node:
                node[p] = {}
            node = node[p]
    return root


def dict_to_ete(node_dict, name="root"):
    t = Tree(name=name)
    for k, v in node_dict.items():
        child = dict_to_ete(v, name=k)
        t.add_child(child)
    return t


tree_dict = build_tree(paths)
ete_tree = dict_to_ete(tree_dict, name="ybyra")  # root

# force leaf for each sample
for path, sample in zip(paths, samples):
    node = ete_tree
    for p in path:
        node = node & p
    node.add_child(Tree(name=sample))

# tree style
ts = TreeStyle()
ts.show_leaf_name = False
ts.branch_vertical_margin = 10

# node styling
for n in ete_tree.traverse():
    nstyle = NodeStyle()
    if not n.is_leaf():
        nstyle["fgcolor"] = "darkred"
        nstyle["size"] = 5
        n.add_face(TextFace(n.name, fsize=10, fgcolor="darkred"), column=0, position="branch-top")
    else:
        nstyle["fgcolor"] = "blue"
        nstyle["size"] = 5
        n.add_face(TextFace(n.name, fsize=10, fgcolor="blue"), column=0, position="branch-right")
    n.set_style(nstyle)

ete_tree.render(output_file, tree_style=ts)

