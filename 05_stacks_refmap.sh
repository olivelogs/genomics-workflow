#!/bin/bash

# 05_stacks_refmap.sh

# Based on my bioinformatics internship work, but in script form

# Runs ref_map.pl on BAM files
# validate_popmap validates popmap file structure 
# Usage: ./05_stacks_refmap.sh <bam_directory> <popmap_file> <output_directory> [threads]

## Arguments ##
INPUT_DIR=$1
POPMAP=$2
OUTPUT_DIR=$3
THREADS=${4:-4}

# check if arguments were provided
if [ -z "$INPUT_DIR" ] || [ -z "$POPMAP" ] || [ -z "$OUTPUT_DIR" ]; then
    echo "Error: Please provide input and output directories and path to popmap"
    echo "Usage: ./05_stacks_refmap.sh <bam_directory> <popmap_file> <output_directory> [threads]"
    exit 1
fi

# Check if input directory exists
if [ ! -d "$INPUT_DIR" ]; then
    echo "Error: Input directory does not exist: $INPUT_DIR"
    exit 1
fi

## Validate popmap files because why not? ##
validate_popmap(){
    local popmap=$1
    local bam_dir=$2
    echo "Validating popmap..."

    # Check if it exists
    if [ ! -f "$popmap" ]; then
        echo "popmap not found: $popmap"
        exit 1
    fi

    # check popmap format 
    while IFS=$'\t' read -r sample pop; do
        # skips blank lines
        [ -z "$sample" ] && continue

        # check if 2 fields filled
        if [ -z "$pop" ]; then
            echo "Error: malformed popmap line (needs tab separation sample and population fields): $sample"
            exit 1
        fi

        # check if corresponding BAM file exists in dir
        if [ ! -f "$bam_dir/${sample}.sorted.bam" ]; then
            echo "Warning: BAM file not found for $sample"
        fi

    done < "$popmap"

    echo "popmap validation complete!"
}

# validate
validate_popmap "$POPMAP" "$INPUT_DIR"


## stacks_refmap ##

# create output dir if it does not exist
mkdir -p "$OUTPUT_DIR"

# run ref_map.pl
ref_map.pl -T "$THREADS" --popmap "$POPMAP" -o "$OUTPUT_DIR" --samples "$INPUT_DIR"

# check if ref_map.pl succeeded
if [ $? -ne 0 ]; then
    echo "error: ref_map.pl failed"
    exit 1
fi

echo "ref_map.pl completed successfully!"