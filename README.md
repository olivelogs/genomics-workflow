# genomics-workflow

This is a re-creation of the workflow used at my 2022 bioinformatics internship!  

**01_fastqc.sh**: runs fastQC on fastq files  
**02_process_radtags.sh**: runs Stacks process_radtags on raw data
*Additional QC steps and diagnostic steps (kmer_filter, stacks_dist_extract) were explored during development but are not included in the final workflow*  
**03_bwa_mem.sh**: runs BWA mem on fq files. Aligns to an existing reference.  
**04_samtools.sh**: runs SAMtools view, sort, index, quickcheck, and flagstat on .sam files from bwa mem  
**05_stacks_refmap.sh**: runs Stacks ref_map.pl on aligned BAM files (gstacks and populations)
**06_vcftools.sh**: uses VCFtools to filter vcf files generated with Stacks. All flags are hard-coded to reflect my original work. Flags are optimized for population genetics analyses with large, high-diversity datasets of distributed populations, preserving rare alleles associated with location-specific adaptations.  
