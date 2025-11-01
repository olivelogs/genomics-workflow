#!/bin/bash

# 05_vcftools.sh

# Based on my bioinformatics internship work, but in script form

# Filters vcf file to select high quality biallelic SNPs, present in at least 65% of samples, with MAF at least 1%, and outputs a new vcf file
    # Why? working with high-diversity wild populations; wanted to keep diversity in the filtering
    # lower MAF preserved rare alleles associated in local adaptations
    # biallelic for compatibility in downstream population genetics analyses
    # all flags are hard-coded to reflect my previous work

# Usage: ./06_vcftools.sh <input_vcf_file_path> <output_prefix>

## Arguments ##
VCF_PATH=$1
OUT_PREFIX=$2

## checks ##
# check if arguments were provided
if [ -z "$VCF_PATH" ] || [ -z "$OUT_PREFIX" ]; then
    echo "Error: Please provide input vcf file and output file prefix"
    echo "Usage: ./06_vcftools.sh <input_vcf_file_path> <output_prefix>"
    exit 1
fi

# Check if input file exists 
if [ ! -f "$VCF_PATH" ]; then
    echo "File not found: $VCF_PATH"
    exit 1
fi

## run VCF filtering ##
echo "Running vcftools on $VCF_PATH..."
vcftools --vcf "$VCF_PATH" --out "$OUT_PREFIX" --maf 0.01 --min-alleles 2 --max-alleles 2 --max-missing 0.65 --remove-indels --recode

# check if VCFtools succeeded
if [ $? -eq 0 ]; then
    echo "vcf file filtered successfully!"
else
    echo "Error: vcf filtering failed"
    exit 1
fi