# genomics-workflow

This is a re-creation of the workflow used at my 2022 bioinformatics internship!  
Used to process bulk RAD-seq data for 176 samples.  

## Scripts  

**01_fastqc.sh**: runs fastQC on fastq files  
**02_process_radtags.sh**: runs Stacks process_radtags on raw data
*Additional QC steps and diagnostic steps (kmer_filter, stacks_dist_extract) were explored during development but are not included in the final workflow*  
**03_bwa_mem.sh**: runs BWA mem on fq files. Aligns to an existing reference.  
**04_samtools.sh**: runs SAMtools view, sort, index, quickcheck, and flagstat on .sam files from bwa mem  
**05_stacks_refmap.sh**: runs Stacks ref_map.pl on aligned BAM files (gstacks and populations)  
**06_vcftools.sh**: uses VCFtools to filter vcf files generated with Stacks. All flags are hard-coded to reflect my original work. Flags are optimized for population genetics analyses with large, high-diversity datasets of distributed populations, preserving rare alleles associated with location-specific adaptations.  

## Next Steps  

**Snakemake**  
i forgot to branch before starting snakemake. lessons were learned. please look at /scripts/ for what is actually functional right now.  
I'm working on turning all of this into a reproducible workflow with Snakemake. Which requires refactoring everything in the scripts, as they loop through directories. So... stay tuned?  
My test data is hidden, you'll need to copy config_template.yaml to config.yaml and fill in your sample info.  

**Landscape Genomics**  
All of the R stuff is coming soon too!  
