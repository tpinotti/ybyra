import sys
import argparse
import os

def parse_yplace(file, min_tree_score):
    best_placement = None
    best_score = float('-inf')
    best_path = ""
    best_ancestral = 0
    ties = []
    all_nodes = {}

    with open(file, 'r') as f:
        next(f)  # Skip header
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 5:
                continue

            node, derived, ancestral, score, path = parts
            derived, ancestral, score = int(derived), int(ancestral), int(score)

            all_nodes[node] = (derived, ancestral, score, path)

            if score >= min_tree_score:
                if score > best_score:
                    best_score = score
                    best_placement = node
                    best_path = path
                    best_ancestral = ancestral
                    ties = [node]
                elif score == best_score:
                    ties.append(node)

    return best_placement, best_score, best_path, best_ancestral, ties, all_nodes


def find_common_parent(ties, all_nodes):
    """Finds the MRCA  with a derived hit for tied nodes."""
    if not ties:
        return None, None

    shortest_path = min(ties, key=lambda x: len(all_nodes[x][3].split('<')))
    paths = [set(all_nodes[node][3].split('<')) for node in ties]
    common_ancestors = set.intersection(*paths)
    valid_ancestors = [a for a in common_ancestors if a in all_nodes]
    sorted_ancestors = sorted(valid_ancestors, key=lambda x: -len(all_nodes[x][3].split('<')))

    for ancestor in sorted_ancestors:
        derived, _, _, _ = all_nodes[ancestor]
        if derived > 0:
            return shortest_path, ancestor

    return shortest_path, None


def has_upstream_support(node, all_nodes, step_size):
    """STEP RULE: check if upstream nodes within the step size have derived hits"""
    if step_size == 0:
        return True  # rule disabled
    if node not in all_nodes:
        return False

    derived, ancestral, score, path = all_nodes[node]
    path_nodes = path.split('<')

    try:
        idx = path_nodes.index(node)
    except ValueError:
        return False

    upstream_nodes = path_nodes[idx + 1 : idx + 1 + step_size]

    for anc in upstream_nodes:
        if anc in all_nodes:
            anc_derived = all_nodes[anc][0]
            if anc_derived > 0:
                return True
    return False


def summarize_sample(file, step_size=5, min_tree_score=10, low_tree_score=50):
    """
    Get the placement of a sample from its yplace file. Returns a dict with the placement,
    or `placement` None and the reason in `fail_flag` if the sample fails, as well as the
    score ties, step rule nopass nodes, and unstable downstream counts of the sample.
    """
    result = {
        "placement": None, "score": None, "flag": None, "path": None, "fail_flag": None,
        "ties": [], "ties_summary": None, "step_rule_nopass": [], "unstable": None
    }
    best_placement, best_score, best_path, best_ancestral, ties, all_nodes = parse_yplace(file, min_tree_score)

    if not best_placement or best_score < min_tree_score:
        result["fail_flag"] = "below_min_tree_score"
        return result

    flag_parts = []

    if best_score < low_tree_score:
        flag_parts.append("low_tree_score")

    if len(ties) > 1:
        flag_parts.append("score_tie")
        result["ties"] = [(node, *all_nodes[node]) for node in ties]
        shortest_path, common_parent = find_common_parent(ties, all_nodes)
        result["ties_summary"] = (shortest_path, common_parent)

        # Always use MRCA for tie resolution
        if common_parent:
            best_placement = common_parent
            best_path = all_nodes[common_parent][3]
            flag_parts.append("most_recent_common_parent")

    step_rule_applied = False

    candidates = sorted(all_nodes.items(), key=lambda x: -x[1][2])

    start_index = next((i for i, (n, _) in enumerate(candidates) if n == best_placement), 0)

    passed = False
    for i in range(start_index, len(candidates)):
        node, (derived, ancestral, score, path) = candidates[i]
        if score < min_tree_score:
            break  # label as fail

        if has_upstream_support(node, all_nodes, step_size):
            if node != best_placement:
                step_rule_applied = True
            best_placement = node
            best_score = score
            best_path = path
            passed = True
            break
        else:
            result["step_rule_nopass"].append((node, score, path))

    if not passed:
        result["fail_flag"] = "below_min_tree_score_after_step_rule"
        return result

    if step_rule_applied:
        flag_parts.append("step_rule")

    derived, ancestral, _, _ = all_nodes[best_placement]
    if ancestral != 0:
        flag_parts.append("unstable_downstream")
        result["unstable"] = (derived, ancestral)

    result["placement"] = best_placement
    result["score"] = best_score
    result["path"] = best_path
    result["flag"] = ";".join(flag_parts) if flag_parts else "..."
    return result


def main(files, step_size=5, min_tree_score=10, low_tree_score=50):
    with open("aggregate.yplace", 'w') as agg, \
         open("score_ties.yplace", 'w') as ties, \
         open("score_ties_summary.yplace", 'w') as ties_summary, \
         open("unstable_downstream.yplace", 'w') as unstable, \
         open("step_rule_nopass.yplace", 'w') as nopass, \
         open("fail.yplace", 'w') as fail:

        agg.write("individual\toptplacement\ttree_score\tflag\ttree_path\n")
        ties.write("individual\tid\tderived\tancestral\ttree_score\ttree_path\n")
        ties_summary.write("individual\tshortest_path_to_root\tmost_recent_common_parent\n")
        unstable.write("individual\tid\tderived\tancestral\ttree_score\ttree_path\n")
        nopass.write("individual\tid\tscore\ttree_path\n")
        fail.write("individual\tflag\n")

        for file in files:
            individual = os.path.basename(file).replace(".yplace", "")

            if not os.path.getsize(file):
                continue

            r = summarize_sample(file, step_size, min_tree_score, low_tree_score)

            for node, derived, ancestral, score, path in r["ties"]:
                ties.write(f"{individual}\t{node}\t{derived}\t{ancestral}\t{score}\t{path}\n")
            if r["ties_summary"]:
                shortest_path, common_parent = r["ties_summary"]
                ties_summary.write(f"{individual}\t{shortest_path}\t{common_parent if common_parent else 'None'}\n")
            for node, score, path in r["step_rule_nopass"]:
                nopass.write(f"{individual}\t{node}\t{score}\t{path}\n")

            if r["placement"] is None:
                fail.write(f"{individual}\t{r['fail_flag']}\n")
                continue

            if r["unstable"]:
                derived, ancestral = r["unstable"]
                unstable.write(f"{individual}\t{r['placement']}\t{derived}\t{ancestral}\t{r['score']}\t{r['path']}\n")
            agg.write(f"{individual}\t{r['placement']}\t{r['score']}\t{r['flag']}\t{r['path']}\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Get optimal placement for multiple yplace files.")
    parser.add_argument("files", nargs='+', help="List of yplace output files to aggregate")
    parser.add_argument(
        "--step-size",
        type=int,
        default=5,
        help="Step size to check upstream nodes; 0 disables the step rule (default: 5)"
    )
    parser.add_argument(
        "--min-tree-score",
        type=int,
        default=10,
        help="Minimum tree score for a placement; samples below fail (default: 10)"
    )
    parser.add_argument(
        "--low-tree-score",
        type=int,
        default=50,
        help="Placements below this tree score are flagged as low_tree_score (default: 50)"
    )
    args = parser.parse_args()
    if args.step_size < 0:
        parser.error("--step-size must be >= 0")

    main(
        args.files,
        step_size=args.step_size,
        min_tree_score=args.min_tree_score,
        low_tree_score=args.low_tree_score
    )
