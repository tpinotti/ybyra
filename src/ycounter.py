import sys
import collections

def count_states(calls_file):
    counts = collections.defaultdict(lambda: {'derived': 0, 'ancestral': 0})
    
    with open(calls_file, 'r') as f:
        next(f)  # Skip header
        for line in f:
            parts = line.strip().split('\t')
            node, state = parts[3], parts[7]
            if state == "derived":
                counts[node]['derived'] += 1
            elif state == "ancestral":
                counts[node]['ancestral'] += 1
    
    return counts

def main(calls_file):
    counts = count_states(calls_file)
    print("id\tderived\tancestral")
    for node, values in counts.items():
        print(f"{node}\t{values['derived']}\t{values['ancestral']}")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python ycounter.py <calls_file>")
        sys.exit(1)
    
    calls_file = sys.argv[1]
    main(calls_file)
