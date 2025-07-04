#!/usr/bin/env python3
"""
FamilyTreeDNA Y-tree JSON to TSV Converter
This script converts FamilyTreeDNA JSON format to the specified TSV format
"""

import json
import csv
import sys
from pathlib import Path


def convert_ftdna_json_to_tsv(json_data):
    """
    Convert FamilyTreeDNA JSON data to TSV format
    
    Args:
        json_data: Either a JSON string or a dictionary containing the FTDNA data
    
    Returns:
        str: TSV formatted string with header and data rows
    """
    # Parse JSON if it's a string
    if isinstance(json_data, str):
        data = json.loads(json_data)
    else:
        data = json_data
    
    # Create maps for efficient lookup
    id_to_name = {}
    id_to_parent_id = {}
    
    # First pass: collect all node names and parent relationships
    for node in data['allNodes'].values():
        id_to_name[node['haplogroupId']] = node['name']
        id_to_parent_id[node['haplogroupId']] = node.get('parentId')
    
    def get_parent_name(parent_id):
        """Get parent haplogroup name from parent ID"""
        if not parent_id:
            return ''
        return id_to_name.get(parent_id, f'ID-{parent_id}')
    
    # Prepare output data
    output_rows = []
    header = ['snpId', 'otherId', 'position', 'mutation', 'anc', 'der', 'id', 'parent']
    output_rows.append(header)
    
    # Process each node
    for node in data['allNodes'].values():
        node_name = node['name']
        parent_name = get_parent_name(node.get('parentId'))
        
        # Process each variant in the node
        variants = node.get('variants', [])
        for variant in variants:
            snp_id = variant.get('variant', '')
            other_id = variant.get('variant', '')
            position = variant.get('position', '')
            ancestral = variant.get('ancestral', '')
            derived = variant.get('derived', '')
            mutation = f"{ancestral}>{derived}" if ancestral and derived else ''
            
            row = [
                snp_id,      # snpId (using variant name)
                other_id,    # otherId (same as snpId)
                str(position) if position else '',  # position
                mutation,    # mutation (anc>der format)
                ancestral,   # anc
                derived,     # der
                node_name,   # id (haplogroup name)
                parent_name  # parent (parent haplogroup name)
            ]
            
            output_rows.append(row)
    
    # Convert to TSV string
    tsv_lines = []
    for row in output_rows:
        tsv_lines.append('\t'.join(row))
    
    return '\n'.join(tsv_lines)


def convert_file(input_file, output_file=None):
    """
    Convert a JSON file to TSV format
    
    Args:
        input_file: Path to input JSON file
        output_file: Path to output TSV file (optional, defaults to input_file with .tsv extension)
    """
    input_path = Path(input_file)
    
    # Set default output file if not provided
    if output_file is None:
        output_file = input_path.with_suffix('.tsv')
    
    # Read input JSON file
    try:
        with open(input_path, 'r', encoding='utf-8') as f:
            json_data = f.read()
    except FileNotFoundError:
        print(f"Error: Input file '{input_file}' not found.")
        return False
    except Exception as e:
        print(f"Error reading input file: {e}")
        return False
    
    # Convert to TSV
    try:
        tsv_output = convert_ftdna_json_to_tsv(json_data)
    except json.JSONDecodeError as e:
        print(f"Error parsing JSON: {e}")
        return False
    except Exception as e:
        print(f"Error during conversion: {e}")
        return False
    
    # Write output TSV file
    try:
        with open(output_file, 'w', encoding='utf-8', newline='') as f:
            f.write(tsv_output)
        print(f"Successfully converted '{input_file}' to '{output_file}'")
        return True
    except Exception as e:
        print(f"Error writing output file: {e}")
        return False


def main():
    """Main function for command line usage"""
    if len(sys.argv) < 2:
        print("Usage: python ftdna_converter.py <input_json_file> [output_tsv_file]")
        print("Example: python ftdna_converter.py ftdna_ytree.json output.tsv")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2] if len(sys.argv) > 2 else None
    
    success = convert_file(input_file, output_file)
    sys.exit(0 if success else 1)


if __name__ == "__main__":
    main()


# Example usage in Python script:
"""
# Method 1: Convert from file
convert_file('ftdna_ytree.json', 'output.tsv')

# Method 2: Convert from JSON string/dict
with open('ftdna_ytree.json', 'r') as f:
    json_data = f.read()

tsv_output = convert_ftdna_json_to_tsv(json_data)
print(tsv_output)

# Method 3: Convert from dictionary
import json
with open('ftdna_ytree.json', 'r') as f:
    data_dict = json.load(f)

tsv_output = convert_ftdna_json_to_tsv(data_dict)
with open('output.tsv', 'w') as f:
    f.write(tsv_output)
"""