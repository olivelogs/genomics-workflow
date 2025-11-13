#!/bin/bash

# 03_bwa_mem.sh

# Based on my bioinformatics internship work, but in script form 

# Runs bwa mem on all fq files in a directory. 
# Usage: ./03_bwa_mem.sh <input directory> <output directory> <path_to_reference> [threads (default 4)]

INPUT_DIR=$1
OUTPUT_DIR=$2
REF_PATH=$3
THREADS=${4:-4} # Defaults to 4 if not provided

# Check if arguments were provided
if [ -z "$INPUT_DIR" ] || [ -z "$OUTPUT_DIR" ] || [ -z "$REF_PATH" ]; then
    echo "Error: Please provide input and output directories and path to reference genome"
    echo "Usage: ./03_bwa_mem.sh <input_directory> <output_directory> <path_to_reference> [threads (default 4)]"
    exit 1
fi

# Check if input directory exists
if [ ! -d "$INPUT_DIR" ]; then
    echo "Error: Input directory does not exist: $INPUT_DIR"
    exit 1
fi

# Check if input directory has files in it
if [ -z "$(ls -A "$INPUT_DIR"/*.fq 2>/dev/null)" ]; then
    echo "Error: no .fq files found in $INPUT_DIR"
    exit 1
fi

# Check if reference.fa exists and has been indexed
if [ ! -f "$REF_PATH" ]; then
    echo "Error: Reference genome not found"
    exit 1
fi

if [ ! -f "${REF_PATH}.bwt" ]; then
    echo "Error: Reference genome doesn't appear to be indexed. Run 'bwa index $REF_PATH' first."
    exit 1
fi

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Run bwa mem on all data files
echo "Running bwa mem on files in $INPUT_DIR with $THREADS threads..." 
for file in "$INPUT_DIR"/*.fq; do
    basename=$(basename "$file" .fq)
    echo "Processing $basename..."
    bwa mem -t "$THREADS" "$REF_PATH" "$file" > "$OUTPUT_DIR/${basename}.sam"

    if [ $? -ne 0 ]; then
        echo "Error: bwa mem failed on $basename"
        exit 1
    fi
done

echo "bwa completed successfully for all files!"