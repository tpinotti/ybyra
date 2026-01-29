import sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

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


# Layout computation: assign y-positions to leaves and compute x from depth
class TreeNode:
    def __init__(self, name, depth=0):
        self.name = name
        self.depth = depth
        self.children = []
        self.y = None
        self.is_sample = False

def dict_to_tree(node_dict, name="root", depth=0):
    t = TreeNode(name, depth)
    for k, v in node_dict.items():
        child = dict_to_tree(v, name=k, depth=depth + 1)
        t.children.append(child)
    return t


tree_dict = build_tree(paths)
root = dict_to_tree(tree_dict, name="ybyra", depth=0)

# attach samples as leaves
for path, sample in zip(paths, samples):
    node = root
    for p in path:
        for c in node.children:
            if c.name == p:
                node = c
                break
    leaf = TreeNode(sample, node.depth + 1)
    leaf.is_sample = True
    node.children.append(leaf)


def assign_y(node, counter=[0]):
    if not node.children:
        node.y = counter[0]
        counter[0] += 1
    else:
        for c in node.children:
            assign_y(c, counter)
        node.y = sum(c.y for c in node.children) / len(node.children)


assign_y(root)


def get_max_depth(node):
    if not node.children:
        return node.depth
    return max(get_max_depth(c) for c in node.children)


def count_leaves(node):
    if not node.children:
        return 1
    return sum(count_leaves(c) for c in node.children)


max_depth = get_max_depth(root)
n_leaves = count_leaves(root)

fig_height = max(6, n_leaves * 0.25)
fig_width = max(10, max_depth * 1.5)
fig, ax = plt.subplots(1, 1, figsize=(fig_width, fig_height))


def draw_tree(node, ax):
    for child in node.children:
        # horizontal line from parent to child's x at parent's y level
        ax.plot([node.depth, child.depth], [child.y, child.y], color="black", lw=0.8)
        # vertical line at parent's x connecting children
        draw_tree(child, ax)

    if len(node.children) > 1:
        ys = [c.y for c in node.children]
        ax.plot([node.depth, node.depth], [min(ys), max(ys)], color="black", lw=0.8)

    if node.is_sample:
        ax.text(node.depth + 0.1, node.y, node.name, fontsize=7,
                va="center", ha="left", color="blue")
        ax.plot(node.depth, node.y, "o", color="blue", markersize=3)
    elif not node.children:
        ax.text(node.depth + 0.1, node.y, node.name, fontsize=7,
                va="center", ha="left", color="darkred")
        ax.plot(node.depth, node.y, "o", color="darkred", markersize=3)
    else:
        ax.text(node.depth, node.y + 0.3, node.name, fontsize=7,
                va="bottom", ha="center", color="darkred")
        ax.plot(node.depth, node.y, "o", color="darkred", markersize=3)


draw_tree(root, ax)
ax.set_xlim(-0.5, max_depth + 3)
ax.set_ylim(-1, n_leaves)
ax.axis("off")
plt.tight_layout()
plt.savefig(output_file, dpi=150, bbox_inches="tight")
plt.close()
