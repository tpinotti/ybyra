#!/bin/bash
# Script to create samples.tsv file from BAM/CRAM files in a directory
# Usage: ./create_samples_tsv.sh PATH OUTPUT_FILENAME

# Check if correct number of arguments provided
if [ $# -ne 2 ]; then
    echo "Usage: $0 PATH OUTPUT_FILENAME"
    echo "Example: $0 /path/to/bam/files samples.tsv"
    exit 1
fi

PATH_TO_BAMS="$1"
OUTPUT_FILE="$2"

# Check if the path exists
if [ ! -d "$PATH_TO_BAMS" ]; then
    echo "Error: Directory '$PATH_TO_BAMS' does not exist"
    exit 1
fi

# Create header
echo -e "sampleId\tbam" > "$OUTPUT_FILE"

# Find all .bam and .cram files in the specified directory only (not subdirectories)
find "$PATH_TO_BAMS" -maxdepth 1 \( -name "*.bam" -o -name "*.cram" \) -type f | while read -r bamfile; do
    # Get the basename without path
    basename_file=$(basename "$bamfile")
    
    # Extract sample ID by removing extensions and common suffixes
    # This handles patterns like: sample.merged.bam, sample.sorted.bam, sample_SE.mapped.bam, etc.
    sample_id=$(echo "$basename_file" | sed 's/\.\(bam\|cram\)$//' | sed 's/\.merged$//' | sed 's/\.sorted$//' | sed 's/_SE\.mapped$//' | sed 's/\.mapped$//')
    
    # Get absolute path to BAM/CRAM file
    absolute_path=$(readlink -f "$bamfile")
    
    # Add to TSV file
    echo -e "${sample_id}\t${absolute_path}"
done >> "$OUTPUT_FILE"

# Count how many samples were found
sample_count=$(tail -n +2 "$OUTPUT_FILE" | wc -l)
echo "Created $OUTPUT_FILE with $sample_count samples"
echo "File contents:"
echo "=============="
cat "$OUTPUT_FILE"