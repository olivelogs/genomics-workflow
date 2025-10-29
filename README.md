# genomics-workflow

This is a re-creation of the workflow used at my 2022 bioinformatics internship!  

**01_fastqc.sh** runs fastQC on fastq files  
**02_process_radtags.sh** runs Stacks process_radtags on raw data
*Additional QC steps and diagnostic steps (kmer_filter, stacks_dist_extract) were explored during development but are not included in the final workflow*  
**03_bwa_mem.sh** runs bwa mem on fq files in a directory  
