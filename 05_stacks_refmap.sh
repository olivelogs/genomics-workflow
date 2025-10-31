#!/bin/bash

# 04_samtools.sh

# Based on my bioinformatics internship work, but in script form

# Runs ref_map.pl on BAM files
# validate_popmap validates popmap file structure 
# Usage: ./05_stacks_refmap.sh <bam_directory> <popmap_file> <output_directory> [threads]

## Arguments ##
INPUT_DIR=$1
POPMAP=$2
OUTPUT_DIR=$3
THREADS=${4:-4}

## Validate popmap files because why not? ##
validate_popmap(){
    local popmap=$1
    local bam_dir=$1
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
    done
}

