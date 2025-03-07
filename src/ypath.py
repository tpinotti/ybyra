import sys

def load_tree(tree_file):
    tree = {}
    with open(tree_file, 'r') as f:
        next(f)  # Skip header
        for line in f:
            child, parent = line.strip().split(',')
            tree[child] = parent
    return tree

def get_path_to_root(tree, node):
    path = []
    while node in tree:
        path.append(node)
        node = tree[node]
        if node == "ybyra":  # Stop at root
            path.append(node)
            break
    return path

def main(tree_file, node, long_format=False):
    tree = load_tree(tree_file)
    if node not in tree:
        print(f"Node {node} not found in tree.")
        return
    
    path = get_path_to_root(tree, node)
    
    if long_format:
        print(" < ".join(path))
    else:
        print("\n".join(path))

if __name__ == "__main__":
    if len(sys.argv) < 3 or len(sys.argv) > 4:
        print("Usage: python ypath.py [--long] <tree_file> <node>")
        sys.exit(1)
    
    long_format = "--long" in sys.argv
    args = [arg for arg in sys.argv[1:] if arg != "--long"]
    
    tree_file = args[0]
    node = args[1]
    
    main(tree_file, node, long_format)

