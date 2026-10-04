"""
Explore how the placements change with the step size of the step rule.

Places all samples at step sizes 1 to a maximum and without the step rule ("off"), using
the same placement as ysummary.py, and compares each placement to a baseline: the placement
at the configured step size, or without the step rule. Writes a table and two plots:
an overview of the changes per step size and by sample depth, and the per-sample placements
for all samples whose placement changes.
"""

import argparse
import csv
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.backends.backend_pdf import PdfPages  # noqa: E402
from matplotlib.patches import Patch, Rectangle  # noqa: E402
import numpy as np  # noqa: E402

import ysummary  # noqa: E402

CLASSES = {
    "identical":                "#4daf4a",
    "upstream 1-5 nodes":       "#c6dbef",
    "upstream 6-10 nodes":      "#6baed6",
    "upstream >10 nodes":       "#08519c",
    "downstream":               "#984ea3",
    "other lineage":            "#e41a1c",
    "failed (low score)":       "#bdbdbd",
    "failed (step rule)":       "#525252",
    "placed (baseline failed)": "#ff7f00",
}
FAIL_CLASSES = {
    "below_min_tree_score": "failed (low score)",
    "below_min_tree_score_after_step_rule": "failed (step rule)",
}
DEPTH_BINS = [0, 0.1, 0.5, 1, 5, float("inf")]
DEPTH_LABELS = ["< 0.1x", "0.1-0.5x", "0.5-1x", "1-5x", "> 5x"]
ROWS_PER_PAGE = 60

# =================================================================================================
#     Input
# =================================================================================================

def read_tree(tree_file):
    with open(tree_file) as f:
        return {row["id"]: row["parent"] for row in csv.DictReader(f)}

def read_snp_info(info_file):
    """Number of SNPs per node, and number of unique SNP positions."""
    snps = {}
    positions = set()
    with open(info_file) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            snps[row["id"]] = snps.get(row["id"], 0) + 1
            positions.add(row["position"])
    return snps, len(positions)

def sample_depth(calls_file, n_positions):
    """Average reads per panel SNP position, from the calls and the calls that did not pass."""
    reads = {}
    nopass_file = calls_file[:-len(".calls")] + ".nopass"
    for path in (calls_file, nopass_file):
        if path == nopass_file and not os.path.exists(path):
            continue
        with open(path) as f:
            for row in csv.DictReader(f, delimiter="\t"):
                # Positions listed for several SNPs share their reads, so count them once.
                reads[row["position"]] = sum(int(x) for x in row["dp4"].split(","))
    return sum(reads.values()) / n_positions

# =================================================================================================
#     Classification
# =================================================================================================

def path_to_root(parent, node):
    path = [node]
    while path[-1] in parent:
        path.append(parent[path[-1]])
    return path

def distance_up(parent, snps, node, ancestor):
    """Nodes and SNPs from `node` up to its `ancestor`."""
    nodes = n_snps = 0
    while node != ancestor:
        n_snps += snps.get(node, 0)
        nodes += 1
        node = parent[node]
    return nodes, n_snps

def classify(parent, snps, result, baseline):
    """Class of a placement relative to the baseline placement, and nodes and SNPs moved up."""
    if result["placement"] is None:
        return FAIL_CLASSES[result["fail_flag"]], None, None
    if baseline["placement"] is None:
        return "placed (baseline failed)", None, None
    p, b = result["placement"], baseline["placement"]
    if p == b:
        return "identical", 0, 0
    if p in path_to_root(parent, b):
        nodes, n_snps = distance_up(parent, snps, b, p)
        if nodes <= 5:
            return "upstream 1-5 nodes", nodes, n_snps
        return ("upstream 6-10 nodes" if nodes <= 10 else "upstream >10 nodes"), nodes, n_snps
    if b in path_to_root(parent, p):
        nodes, n_snps = distance_up(parent, snps, p, b)
        return "downstream", -nodes, -n_snps
    return "other lineage", None, None

def outcome(result):
    """What is compared to tell whether a placement changed: the node, or the reason for failing."""
    return result["placement"] or result["fail_flag"]

# =================================================================================================
#     Plots
# =================================================================================================

def step_label(step):
    return "off" if step == 0 else str(step)

def step_positions(steps):
    """x positions of the steps, with "off" (step 0) on the right, after a gap."""
    return [float(s) if s > 0 else max(steps) + 1.5 for s in steps]

def plot_overview(rows, steps, baseline_step, baseline_name, out_file):
    xpos = step_positions(steps)
    n_samples = len({r["sample"] for r in rows})
    fig, (ax_bars, ax_depth) = plt.subplots(1, 2, figsize=(12, 4.8))

    # Fraction of samples per class and step.
    bottom = np.zeros(len(steps))
    for cls, color in CLASSES.items():
        frac = np.array([
            sum(1 for r in rows if r["step"] == s and r["class"] == cls) / n_samples for s in steps
        ])
        if frac.any():
            ax_bars.bar(xpos, frac, bottom=bottom, color=color, width=0.8, label=cls)
        bottom += frac
    ax_bars.set_ylim(0, 1)
    ax_bars.set_ylabel(f"fraction of samples (n = {n_samples})")
    ax_bars.legend(frameon=False, fontsize=8, loc="upper left", bbox_to_anchor=(1.0, 1.0))

    # Fraction of changed placements per step, by sample depth, as a gradient over the depth bins.
    depth_bins = np.digitize([r["depth"] for r in rows], DEPTH_BINS[1:-1])
    colors = plt.get_cmap("viridis")(np.linspace(0, 0.85, len(DEPTH_LABELS)))
    for b, (label, color) in enumerate(zip(DEPTH_LABELS, colors)):
        sub = [r for r, rb in zip(rows, depth_bins) if rb == b]
        if not sub:
            continue
        n = len({r["sample"] for r in sub})
        frac = [
            np.mean([r["changed"] for r in sub if r["step"] == s]) for s in steps
        ]
        ax_depth.plot(xpos[:-1], frac[:-1], marker="o", color=color, label=f"{label} (n = {n})")
        ax_depth.plot(xpos[-1:], frac[-1:], marker="o", color=color)
    ax_depth.set_ylim(-0.02, 1.02)
    ax_depth.set_ylabel("fraction of placements changed")
    ax_depth.legend(title="average reads per panel SNP", frameon=False, fontsize=8, title_fontsize=8)

    for ax in (ax_bars, ax_depth):
        ax.set_xticks(xpos, [step_label(s) for s in steps])
        ax.set_xlabel("step size")
        ax.axvline(xpos[steps.index(baseline_step)], color="#999999", lw=0.8, zorder=0)
    fig.suptitle(f"Placements compared to the baseline ({baseline_name})")
    fig.tight_layout()
    fig.savefig(out_file)
    plt.close(fig)

def plot_samples(rows, steps, baseline_step, baseline_name, out_file):
    xpos = list(range(len(steps)))
    samples = {}
    for r in rows:
        samples.setdefault(r["sample"], {})[r["step"]] = r
    changed = [s for s, by_step in samples.items() if any(r["changed"] for r in by_step.values())]
    changed.sort(key=lambda s: -samples[s][baseline_step]["depth"])
    note = (
        f"{len(changed)} of {len(samples)} samples shown; the other {len(samples) - len(changed)} "
        f"have the same placement at all step sizes. Baseline: {baseline_name}."
    )
    legend = [Patch(color=c, label=cls) for cls, c in CLASSES.items()]

    with PdfPages(out_file) as pdf:
        if not changed:
            fig = plt.figure(figsize=(8, 2))
            fig.text(0.5, 0.5, note, ha="center", va="center", wrap=True)
            pdf.savefig(fig)
            plt.close(fig)
            return

        for start in range(0, len(changed), ROWS_PER_PAGE):
            page = changed[start:start + ROWS_PER_PAGE]
            # Fixed margin at the bottom for legend and note, independent of the number of rows.
            height = 2.2 + 0.25 * len(page)
            fig, ax = plt.subplots(figsize=(2.5 + 1.1 * len(steps), height))
            for y, sample in enumerate(page):
                baseline = samples[sample][baseline_step]
                for x, step in zip(xpos, steps):
                    r = samples[sample][step]
                    ax.add_patch(Rectangle(
                        (x - 0.5, y - 0.5), 1, 1, facecolor=CLASSES[r["class"]],
                        edgecolor="white", lw=0.5
                    ))
                    # Placement names in the baseline column, and where they differ from it.
                    if r["placement"] and (step == baseline_step or r["changed"]):
                        dark = r["class"] in ("upstream >10 nodes", "downstream", "other lineage")
                        ax.text(x, y, r["placement"], ha="center", va="center", fontsize=6,
                                color="white" if dark else "black")
            col = steps.index(baseline_step)
            ax.add_patch(Rectangle(
                (col - 0.5, -0.5), 1, len(page), fill=False, edgecolor="black", lw=1.5
            ))
            ax.set_xlim(-0.5, len(steps) - 0.5)
            ax.set_ylim(len(page) - 0.5, -0.5)
            ax.set_xticks(xpos, [step_label(s) for s in steps])
            ax.xaxis.tick_top()
            ax.set_xlabel("step size")
            ax.xaxis.set_label_position("top")
            ax.set_yticks(range(len(page)), [
                f"{s} ({samples[s][baseline_step]['depth']:.2g}x)" for s in page
            ], fontsize=7)
            ax.tick_params(length=0)
            for spine in ax.spines.values():
                spine.set_visible(False)
            margin = 1.0 / height
            fig.legend(handles=legend, loc="lower center", bbox_to_anchor=(0.5, 0.3 * margin),
                       ncol=5, frameon=False, fontsize=7)
            fig.text(0.5, 0.1 * margin, note, ha="center", va="bottom", fontsize=7)
            fig.tight_layout(rect=(0, margin, 1, 1))
            pdf.savefig(fig)
            plt.close(fig)

# =================================================================================================
#     Main
# =================================================================================================

def main(args):
    parent = read_tree(args.tree)
    snps, n_positions = read_snp_info(args.info)
    max_step = max(args.max_step, args.step_size)
    steps = list(range(1, max_step + 1)) + [0]
    baseline_step = args.step_size if args.baseline == "step_size" else 0
    baseline_name = (
        f"step size {args.step_size}" if args.baseline == "step_size" and args.step_size > 0
        else "no step rule"
    )
    calls = {os.path.basename(f)[:-len(".calls")]: f for f in args.calls}

    rows = []
    for yplace_file in args.yplace:
        sample = os.path.basename(yplace_file).replace(".yplace", "")
        # Same as ysummary.py, which skips empty files.
        if not os.path.getsize(yplace_file):
            continue
        depth = sample_depth(calls[sample], n_positions)
        results = {
            step: ysummary.summarize_sample(
                yplace_file, step, args.min_tree_score, args.low_tree_score
            ) for step in steps
        }
        baseline = results[baseline_step]
        for step, result in results.items():
            cls, nodes, n_snps = classify(parent, snps, result, baseline)
            rows.append({
                "sample": sample, "depth": depth, "step": step,
                "baseline": step == baseline_step,
                "placement": result["placement"], "tree_score": result["score"],
                "flag": result["flag"] or result["fail_flag"],
                "class": cls, "nodes_up": nodes, "snps_up": n_snps,
                "changed": outcome(result) != outcome(baseline),
            })

    os.makedirs(args.out_dir, exist_ok=True)
    columns = [
        "sample", "depth", "step", "baseline", "placement", "tree_score", "flag",
        "class", "nodes_up", "snps_up"
    ]
    with open(os.path.join(args.out_dir, "step_size.tsv"), "w") as out:
        out.write("\t".join(columns) + "\n")
        for r in rows:
            values = dict(r, depth=f"{r['depth']:.4g}", step=step_label(r["step"]))
            out.write("\t".join("" if values[c] is None else str(values[c]) for c in columns) + "\n")

    plot_overview(rows, steps, baseline_step, baseline_name, os.path.join(args.out_dir, "overview.pdf"))
    plot_samples(rows, steps, baseline_step, baseline_name, os.path.join(args.out_dir, "samples.pdf"))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Explore how placements change with the step size.")
    parser.add_argument("--tree", required=True, help="Tree file (CSV with id,parent)")
    parser.add_argument("--info", required=True, help="Haplogroup SNP info table")
    parser.add_argument("--yplace", nargs="+", required=True, help="yplace files of the samples")
    parser.add_argument("--calls", nargs="+", required=True, help="calls files of the samples")
    parser.add_argument("--step-size", type=int, default=5, help="Configured step size (default: 5)")
    parser.add_argument(
        "--baseline", choices=["step_size", "off"], default="step_size",
        help="Compare to the placements at the configured step size, or without the step rule"
    )
    parser.add_argument("--max-step", type=int, default=10, help="Largest step size to explore (default: 10)")
    parser.add_argument("--min-tree-score", type=int, default=10, help="As for ysummary.py (default: 10)")
    parser.add_argument("--low-tree-score", type=int, default=50, help="As for ysummary.py (default: 50)")
    parser.add_argument("--out-dir", default="step_size", help="Output directory (default: step_size)")
    args = parser.parse_args()
    if args.step_size < 0 or args.max_step < 1:
        parser.error("--step-size must be >= 0, and --max-step >= 1")
    main(args)
