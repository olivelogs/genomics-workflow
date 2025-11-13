#!/bin/bash

# 01_fastqc.sh

# Based on my bioinformatics internship work, but in script form 
# Runs FastQC on all FASTQ files in a directory


# Purpose: Run FastQC on all FASTQ files in a directory
# Usage: ./01_fastqc.sh <input_directory> <output_directory> [threads (default 4)]

# Get arguments
INPUT_DIR=$1
OUTPUT_DIR=$2
THREADS=${3:-4} # Defaults to 4

# Check if arguments were provided
if [ -z "$INPUT_DIR" ] || [ -z "$OUTPUT_DIR" ]; then
    echo "Error: Please provide input and output directories"
    echo "Usage: ./01_fastqc.sh <input_directory> <output_directory> [threads (default 4)]"
    exit 1
fi

# Check if input directory exists
if [ ! -d "$INPUT_DIR" ]; then
    echo "Error: Input directory does not exist: $INPUT_DIR"
    exit 1
fi

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Run FastQC on all FASTQ files
echo "Running FastQC on files in $INPUT_DIR..."
fastqc -t "$THREADS" "$INPUT_DIR"/*.fq -o "$OUTPUT_DIR"

# Check if FastQC succeeded
if [ $? -eq 0 ]; then
    echo "FastQC completed successfully!"
else
    echo "Error: FastQC failed"
    exit 1
fi