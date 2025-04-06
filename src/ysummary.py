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

            if score >= 15:
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
    """Finds the most recent common ancestor with a derived hit for tied nodes."""
    if not ties:
        return None, None

    # Find the node with the shortest path to root
    shortest_path = min(ties, key=lambda x: len(all_nodes[x][3].split('<')))

    # Get the full paths of all tied nodes
    paths = [set(all_nodes[node][3].split('<')) for node in ties]

    # Find common ancestors by intersecting paths
    common_ancestors = set.intersection(*paths)

    # Filter to only nodes present in all_nodes to avoid KeyError
    valid_ancestors = [a for a in common_ancestors if a in all_nodes]

    # Sort by path length (most recent first)
    sorted_ancestors = sorted(valid_ancestors, key=lambda x: -len(all_nodes[x][3].split('<')))

    # Find the most recent ancestor with a derived hit
    for ancestor in sorted_ancestors:
        derived, _, _, _ = all_nodes[ancestor]
        if derived > 0:
            return shortest_path, ancestor

    return shortest_path, None  # If no ancestor has a derived hit, return None

def main(files):
    with open("aggregate.yplace", 'w') as out:
        out.write("individual\toptplacement\ttree_score\tflag\ttree_path\n")

    with open("scoreties.yplace", 'w') as out:
        out.write("individual\tid\tderived\tancestral\ttree_score\ttree_path\n")

    with open("scoreties_summary.yplace", 'w') as out:
        out.write("individual\tshortest_path_to_root\tmost_recent_common_parent\n")

    with open("unstabledownstream.yplace", 'w') as out:
        out.write("individual\tid\tderived\tancestral\ttree_score\ttree_path\n")

    for file in files:
        individual = os.path.basename(file).replace(".yplace", "")
        
        if not os.path.getsize(file):
            continue  # Ignore empty files

        best_placement, best_score, best_path, best_ancestral, ties, all_nodes = parse_yplace(file)

        if not best_placement or best_score < 15:
            continue  # Ignore placements with score < 15

        flag_parts = []  

        if len(ties) > 1:
            flag_parts.append("score_tie")
            with open("scoreties.yplace", 'a') as out:
                for node in ties:
                    derived, ancestral, score, path = all_nodes[node]
                    out.write(f"{individual}\t{node}\t{derived}\t{ancestral}\t{score}\t{path}\n")

            shortest_path, common_parent = find_common_parent(ties, all_nodes)

            if shortest_path:
                flag_parts.append("shortest_path_to_root")
            if common_parent and common_parent == shortest_path:
                flag_parts.append("most_recent_common_parent")

            with open("scoreties_summary.yplace", 'a') as out:
                out.write(f"{individual}\t{shortest_path}\t{common_parent if common_parent else 'None'}\n")

            # Update best placement to be the shortest_path_to_root
            best_placement = shortest_path
            best_path = all_nodes[shortest_path][3]  

        # Detect unstable downstream cases: best placement has ancestral ≠ 0
        if best_ancestral != 0:
            flag_parts.append("unstable_downstream")
            with open("unstabledownstream.yplace", 'a') as out:
                out.write(f"{individual}\t{best_placement}\t{all_nodes[best_placement][0]}\t{best_ancestral}\t{best_score}\t{best_path}\n")

        flag = ";".join(flag_parts) if flag_parts else "..."

        with open("aggregate.yplace", 'a') as out:
            out.write(f"{individual}\t{best_placement}\t{best_score}\t{flag}\t{best_path}\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Summarize optimal placements from multiple yplace output files.")
    parser.add_argument("files", nargs='+', help="List of yplace output files to aggregate")
    args = parser.parse_args()
    
    main(args.files)

