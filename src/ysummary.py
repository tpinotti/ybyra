import sys
import argparse
import os

def parse_yplace(file):
    best_placement = None
    best_score = float('-inf')
    best_path = ""
    
    with open(file, 'r') as f:
        next(f)  # Skip header
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 5:
                continue
            
            node, derived, ancestral, score, path = parts
            derived, ancestral, score = int(derived), int(ancestral), int(score)
            
            if ancestral == 0 and score > best_score and score >= 5:  # Ensure score >= 5
                best_score = score
                best_placement = node
                best_path = path
    
    return best_placement, best_score, best_path

def main(files):
    with open("aggregate.yplace", 'w') as out:
        out.write("individual\toptplacement\ttree_score\ttree_path\n")
        
        for file in files:
            individual = os.path.basename(file).replace(".yplace", "")
            best_placement, best_score, best_path = parse_yplace(file)
            
            if best_placement:
                out.write(f"{individual}\t{best_placement}\t{best_score}\t{best_path}\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Summarize optimal placements from multiple yplace output files.")
    parser.add_argument("files", nargs='+', help="List of yplace output files to aggregate")
    args = parser.parse_args()
    
    main(args.files)

