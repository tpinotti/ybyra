#!/usr/bin/env python3
"""
Convert TSV output to BED file and tree file
Creates two files from the FTDNA TSV output:
1. BED file with chromosome positions
2. Tree file with parent-child relationships
"""

import csv
import sys
from pathlib import Path


def create_bed_file(tsv_file, output_bed_file=None):
    """
    Create BED file from TSV with chromosome positions
    
    Args:
        tsv_file: Path to input TSV file
        output_bed_file: Path to output BED file (optional)
    """
    if output_bed_file is None:
        output_bed_file = Path(tsv_file).with_suffix('.bed')
    
    positions = set()
    
    # Read TSV file and collect unique positions
    with open(tsv_file, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f, delimiter='\t')
        
        for row in reader:
            position = row.get('position', '').strip()
            if position and position.isdigit():
                pos_int = int(position)
                # BED format uses 0-based coordinates, so position-1 to position
                positions.add((pos_int - 1, pos_int))
    
    # Sort positions
    sorted_positions = sorted(positions)
    
    # Write BED file
    with open(output_bed_file, 'w', encoding='utf-8') as f:
        # Write header
        f.write("#chrY\t-1\tposition\n")
        
        # Write positions
        for start, end in sorted_positions:
            f.write(f"chrY\t{start}\t{end}\n")
    
    print(f"Created BED file: {output_bed_file} with {len(sorted_positions)} positions")
    return output_bed_file


def create_tree_file(tsv_file, output_tree_file=None, root_name="ybyra"):
    """
    Create tree file from TSV with parent-child relationships
    
    Args:
        tsv_file: Path to input TSV file
        output_tree_file: Path to output tree file (optional)
        root_name: Name for the artificial root (default: "ybyra")
    """
    if output_tree_file is None:
        output_tree_file = Path(tsv_file).with_suffix('.tree')
    
    # Collect all haplogroups and their parents
    haplogroups = set()
    parent_relationships = {}
    
    with open(tsv_file, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f, delimiter='\t')
        
        for row in reader:
            haplogroup_id = row.get('id', '').strip()
            parent = row.get('parent', '').strip()
            
            if haplogroup_id:
                haplogroups.add(haplogroup_id)
                if parent:
                    parent_relationships[haplogroup_id] = parent
                    haplogroups.add(parent)
    
    # Find root haplogroups (those without parents in our dataset)
    roots = []
    for haplogroup in haplogroups:
        if haplogroup not in parent_relationships:
            roots.append(haplogroup)
    
    # Write tree file
    with open(output_tree_file, 'w', encoding='utf-8') as f:
        # Write header
        f.write("id,parent\n")
        
        # Write root haplogroups pointing to artificial root
        for root in sorted(roots):
            f.write(f"{root},{root_name}\n")
        
        # Write all parent-child relationships
        for child, parent in sorted(parent_relationships.items()):
            f.write(f"{child},{parent}\n")
    
    print(f"Created tree file: {output_tree_file} with {len(haplogroups)} haplogroups")
    print(f"Root haplogroups: {sorted(roots)}")
    return output_tree_file


def process_tsv_file(tsv_file, bed_file=None, tree_file=None, root_name="ybyra"):
    """
    Process TSV file to create both BED and tree files
    
    Args:
        tsv_file: Path to input TSV file
        bed_file: Path to output BED file (optional)
        tree_file: Path to output tree file (optional)
        root_name: Name for the artificial root in tree file
    """
    tsv_path = Path(tsv_file)
    
    # Set default output filenames if not provided
    if bed_file is None:
        bed_file = tsv_path.with_name(f"{tsv_path.stem}.bed")
    
    if tree_file is None:
        tree_file = tsv_path.with_name(f"{tsv_path.stem}.tree")
    
    # Create both files
    print(f"Processing TSV file: {tsv_file}")
    bed_output = create_bed_file(tsv_file, bed_file)
    tree_output = create_tree_file(tsv_file, tree_file, root_name)
    
    return bed_output, tree_output


def main():
    """Main function for command line usage"""
    if len(sys.argv) < 2:
        print("Usage: python tsv_to_bed_tree.py <input_tsv_file> [bed_file] [tree_file] [root_name]")
        print("Example: python tsv_to_bed_tree.py output.tsv positions.bed haplotree.tree ybyra")
        print("\nIf bed_file and tree_file are not specified, they will be created with the same")
        print("base name as the input file with .bed and .tree extensions.")
        sys.exit(1)
    
    tsv_file = sys.argv[1]
    bed_file = sys.argv[2] if len(sys.argv) > 2 else None
    tree_file = sys.argv[3] if len(sys.argv) > 3 else None
    root_name = sys.argv[4] if len(sys.argv) > 4 else "ybyra"
    
    try:
        bed_output, tree_output = process_tsv_file(tsv_file, bed_file, tree_file, root_name)
        print(f"\nSuccessfully created:")
        print(f"  BED file: {bed_output}")
        print(f"  Tree file: {tree_output}")
    except FileNotFoundError:
        print(f"Error: Input file '{tsv_file}' not found.")
        sys.exit(1)
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()


# Example usage in Python script:
"""
# Method 1: Create both files with default names
process_tsv_file('output.tsv')

# Method 2: Specify custom output filenames
process_tsv_file('output.tsv', 'positions.bed', 'haplotree.tree')

# Method 3: Create files separately
create_bed_file('output.tsv', 'positions.bed')
create_tree_file('output.tsv', 'haplotree.tree', 'ybyra')
"""