#!/bin/bash

# 02_process_radtags.sh

# Based on my bioinformatics internship work, but in script form 
# To reflect actual work done, the script contains the flags used in the original pipeline 

# Runs process_radtags (Stacks) on all raw reads data in a directory
# Usage: ./02_process_radtags.sh <input_directory> <output_directory> <barcodes.txt optional>


# Get arguments
INPUT_DIR=$1
OUTPUT_DIR=$2
BARCODES=$3 # optional

# Check if arguments were provided
if [ -z "$INPUT_DIR" ] || [ -z "$OUTPUT_DIR" ]; then
    echo "Error: Please provide input and output directories"
    echo "Usage: ./02_process_radtags.sh <input_directory> <output_directory> <barcodes.txt optional>"
    exit 1
fi

# Check if input directory exists
if [ ! -d "$INPUT_DIR" ]; then
    echo "Error: Input directory does not exist: $INPUT_DIR"
    exit 1
fi

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Run process_radtags on all data files
echo "Running process_radtags on files in $INPUT_DIR..."
if [ -n "$BARCODES" ]; then
    process_radtags -p "$INPUT_DIR" -o "$OUTPUT_DIR" -b "$BARCODES" -e ecoRI -c -r -q --disable-rad-check
else
    process_radtags -p "$INPUT_DIR" -o "$OUTPUT_DIR" -e ecoRI -c -r -q --disable-rad-check
fi


# Check if process_radtags succeeded
if [ $? -eq 0 ]; then
    echo "process_radtags completed successfully!"
else
    echo "Error: process_radtags failed"
    exit 1
fi