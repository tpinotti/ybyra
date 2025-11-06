import sys
import argparse
import os

def parse_yplace(file):
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

            if score >= 10:
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


def main(files, step_size=5):
    with open("aggregate.yplace", 'w') as out:
        out.write("individual\toptplacement\ttree_score\tflag\ttree_path\n")

    with open("scoreties.yplace", 'w') as out:
        out.write("individual\tid\tderived\tancestral\ttree_score\ttree_path\n")

    with open("scoreties_summary.yplace", 'w') as out:
        out.write("individual\tshortest_path_to_root\tmost_recent_common_parent\n")

    with open("unstabledownstream.yplace", 'w') as out:
        out.write("individual\tid\tderived\tancestral\ttree_score\ttree_path\n")

    with open(f"step{step_size}nopass.yplace", 'w') as out:
        out.write("individual\tid\tscore\ttree_path\n")

    with open("fail.yplace", 'w') as out:
        out.write("individual\tflag\n")

    for file in files:
        individual = os.path.basename(file).replace(".yplace", "")

        if not os.path.getsize(file):
            continue

        best_placement, best_score, best_path, best_ancestral, ties, all_nodes = parse_yplace(file)

        if not best_placement or best_score < 10:
            with open("fail.yplace", 'a') as out:
                out.write(f"{individual}\tscore_below_10\n")
            continue

        flag_parts = []

        if best_score < 50:
            flag_parts.append("tree_score_below_50")

        if len(ties) > 1:
            flag_parts.append("score_tie")

            with open("scoreties.yplace", 'a') as out:
                for node in ties:
                    derived, ancestral, score, path = all_nodes[node]
                    out.write(f"{individual}\t{node}\t{derived}\t{ancestral}\t{score}\t{path}\n")

            shortest_path, common_parent = find_common_parent(ties, all_nodes)

            with open("scoreties_summary.yplace", 'a') as out:
                out.write(f"{individual}\t{shortest_path}\t{common_parent if common_parent else 'None'}\n")

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
            if score < 10:
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
                with open(f"step{step_size}nopass.yplace", 'a') as out:
                    out.write(f"{individual}\t{node}\t{score}\t{path}\n")

        if not passed:
            with open("fail.yplace", 'a') as out:
                out.write(f"{individual}\tscore_below_10\n")
            continue

        if step_rule_applied:
            flag_parts.append("step_rule")

        derived, ancestral, _, _ = all_nodes[best_placement]
        if ancestral != 0:
            flag_parts.append("unstable_downstream")
            with open("unstabledownstream.yplace", 'a') as out:
                out.write(f"{individual}\t{best_placement}\t{derived}\t{ancestral}\t{best_score}\t{best_path}\n")

        flag = ";".join(flag_parts) if flag_parts else "..."

        with open("aggregate.yplace", 'a') as out:
            out.write(f"{individual}\t{best_placement}\t{best_score}\t{flag}\t{best_path}\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Get optimal placement for multiple yplace files.")
    parser.add_argument("files", nargs='+', help="List of yplace output files to aggregate")
    parser.add_argument(
        "--step-size",
        type=int,
        default=5,
        help="Step size to check upstream nodes (default: 5)"
    )
    args = parser.parse_args()

    main(args.files, step_size=args.step_size)
