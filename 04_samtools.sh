#!/bin/bash

# 04_samtools.sh

# Based on my bioinformatics internship work, but in script form

# Processes SAM files from BWA alignment using samtools
# Usage: ./04_samtools.sh <input_directory> <output_directory> [threads]


INPUT_DIR=$1
OUTPUT_DIR=$2
THREADS=${3:-4} # Defaults to 4 if not provided

# Check if arguments were provided
if [ -z "$INPUT_DIR" ] || [ -z "$OUTPUT_DIR" ]; then
    echo "Error: Please provide input and output directories and path to reference genome"
    echo "Usage: ./04_samtools.sh <input_directory> <output_directory> [threads]"
    exit 1
fi

# Check if input directory exists
if [ ! -d "$INPUT_DIR" ]; then
    echo "Error: Input directory does not exist: $INPUT_DIR"
    exit 1
fi

# Create output directory & stats subdirectory if it doesn't exist
mkdir -p "$OUTPUT_DIR"
mkdir -p "$OUTPUT_DIR/stats"

# Run samtools on .sam files
echo "Processing SAM files in $INPUT_DIR with $THREADS threads..."
for sam_file in "$INPUT_DIR"/*.sam; do
    basename=$(basename "$sam_file" .sam)
    echo "Processing $basename..."
    
    # Convert SAM to BAM (temp)
    samtools view -@ "$THREADS" -b "$sam_file" > "$OUTPUT_DIR/${basename}.bam"
    
    # Sort BAM (keep)
    samtools sort -@ "$THREADS" "$OUTPUT_DIR/${basename}.bam" -o "$OUTPUT_DIR/${basename}.sorted.bam"
    if [ $? -ne 0 ]; then
        echo "Error: sort failed on $basename"
        exit 1
    fi
    
    # Delete unsorted BAM (saves space)
    rm "$OUTPUT_DIR/${basename}.bam"
    
    # Index the sorted BAM
    samtools index "$OUTPUT_DIR/${basename}.sorted.bam"
    
    # Verify integrity
    samtools quickcheck "$OUTPUT_DIR/${basename}.sorted.bam"
    if [ $? -ne 0 ]; then
        echo "Error: BAM file is corrupted for $basename"
        exit 1
    fi
    
    # Get statistics and put them in stats/ subdirectory
    samtools flagstat "$OUTPUT_DIR/${basename}.sorted.bam" > "$OUTPUT_DIR/stats/${basename}.flagstat.txt"
done

echo "All SAM files processed successfully!"