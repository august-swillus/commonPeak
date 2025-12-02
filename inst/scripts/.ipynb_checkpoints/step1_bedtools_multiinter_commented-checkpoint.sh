#!/bin/bash

# bedtools_multiinter.sh
# Usage: bash bedtools_multiinter.sh beds_manifest.tsv output_dir

# Get arguments from R
bed_manifest="$1"      # BED manifest TSV from R
output_dir="$2"        # Output directory from R

# Ensure output directory exists
mkdir -p "$output_dir"

# Prepare list of BED files from manifest
# Assumes manifest has a header: group, sample, path, type
bed_files=($(tail -n +2 "$bed_manifest" | cut -f3)) # take bed_manifest df, starting from row extract col 3, the filepaths, write them into an array bed_files

# Output files
multiinter_output_file="${output_dir}/multiinter_total_output.bed"
final_peaks_file="${output_dir}/peak_in_all_samples.txt"

echo "Running bedtools multiinter on ${#bed_files[@]} files..."

# ------- Run bedtools multiinter --------
bedtools multiinter -i "${bed_files[@]}" > "$multiinter_output_file" # @ extract all elements from array into space separted bash list 

if [ $? -eq 0 ]; then
    echo "bedtools multiinter completed successfully. Output written to: $multiinter_output_file"
    
    # Determine number of BED files (to filter peaks present in all)
    num_files=${#bed_files[@]} # #bed_files[@] extracts number of elemenmts in the array 
    echo "Filtering for peaks present in all $num_files files (column 4 == $num_files)..."
    
    awk -v n="$num_files" '$4 == n' "$multiinter_output_file" > "$final_peaks_file"
    
    if [ $? -eq 0 ]; then
        echo "Filtered peaks (in all $num_files) written to: $final_peaks_file"
    else
        echo "Error during awk filtering."
        exit 1
    fi
else
    echo "Error running bedtools multiinter."
    exit 1
fi

echo "Script finished."
