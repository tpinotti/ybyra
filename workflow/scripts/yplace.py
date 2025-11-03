import sys
import collections
import argparse

def load_tree(tree_file):
    tree = {}
    with open(tree_file, 'r') as f:
        next(f)  # Skip header
        for line in f:
            child, parent = line.strip().split(',')
            tree[child] = parent
    return tree

def count_states(calls_file, no_damage):
    counts = collections.defaultdict(lambda: {'derived': 0, 'ancestral': 0})
    
    with open(calls_file, 'r') as f:
        next(f)  # Skip header
        for line in f:
            parts = line.strip().split('\t')
            node, state, damage = parts[3], parts[7], parts[8]
            
            if no_damage and damage != "no":
                continue  # Skip entries where damage is not "no" if --no-damage is set
            
            if state == "derived":
                counts[node]['derived'] += 1
            elif state == "ancestral":
                counts[node]['ancestral'] += 1
    
    return counts

def get_path_to_root(tree, node):
    path = []
    while node in tree:
        path.append(node)
        node = tree[node]
        if node == "ybyra":  # Stop at root
            path.append(node)
            break
    return path

def calculate_tree_score(tree, counts):
    scores = {}
    
    for node in counts:
        path = get_path_to_root(tree, node)
        score = sum(counts[n]['derived'] - counts[n]['ancestral'] for n in path if n in counts)
        scores[node] = (counts[node]['derived'], counts[node]['ancestral'], score, "<".join(path))
    
    return scores

def main(tree_file, calls_file, output_file, no_damage):
    tree = load_tree(tree_file)
    counts = count_states(calls_file, no_damage)
    scores = calculate_tree_score(tree, counts)
    
    sorted_scores = sorted(scores.items(), key=lambda x: x[1][2], reverse=True)
    
    with open(output_file, 'w') as out:
        out.write("id\tderived\tancestral\ttree_score\ttree_path\n")
        for node, (derived, ancestral, score, path) in sorted_scores:
            if derived > 0:
                out.write(f"{node}\t{derived}\t{ancestral}\t{score}\t{path}\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compute placement scores for individuals in a phylogenetic tree.")
    parser.add_argument("tree_file", help="Path to the tree file (TSV format: id,parent)")
    parser.add_argument("calls_file", help="Path to the calls file (TSV format: position, snpId, mutation, id, parent, coverage, dp4, state, damage)")
    parser.add_argument("output_file", help="Path to save the output file")
    parser.add_argument("--no-damage", action="store_true", help="Exclude SNPs with 'damage' marked as 'yes' from the calculations")
    
    args = parser.parse_args()
    
    main(args.tree_file, args.calls_file, args.output_file, args.no_damage)

