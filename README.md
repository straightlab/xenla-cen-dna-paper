
This repo contains the pipeline used for data analysis in the paper 'Identification and characterization of centromeric sequences in Xenopus laevis' which can be found at https://www.biorxiv.org/content/10.1101/2020.06.23.167643v1 

###Pipeline Overview:

-Takes unprocessed, zipped fastqs as input. Removes PCR duplicates using clumpify.sh and trims adapters using trimmomatic. 
-Processed reads are then used as input for KMC which generates a k-mer database for each fastq. K-mer databases are exported to text files as well as fasta files. 
-The abundance of each k-mer in each k-mer database is then counted and normalized to the number of basepairs in the processed fastq file. 
-An enrichment ratio is then calculated for each k-mer by dividing the number of times it is found in dataset1 by the number of times it is found in dataset2 (where dataset1 and dataset2 are specified in the config file as 'pairing'). 
-K-mers are then filtered based on enrichment value cutoff (specified in config). Cutoffs are based on the number of median absolute deviations from the median enrichment value for the given pairing. K-mers with enrichment values above the specified cutoff are written to a fasta file.
-Filtered k-mers are then used to select reads from a specified fastq (in config file: 'your_fav_seqdata'). Reads containing at least one k-mer from input k-mer file are outputted to a fasta file. This file is then converted into a k-mer table file where each row is a read and each column is a k-mer ID number. Each line contains either a 1 or 0 denoting the presence or absence of the corresponding k-mer in that read.


###How to run it:

-edit config file:
--'WORKING_DIR' should be the directory where you will run the pipeline. All output directories and files will be within this directory. The pipeline file as well as the config file should also be in this directory.
--'samples:' specifies the entire filepath for raw fastqs. The first indentation level indicates the name of the parent direcrtory that will be created. The second indicates the {basename}. The third idicates the readnum and filepath. For single end reads, each readnum is r1.
--'your_fav_seqdata:' specifies the fastqs that will be used as input in BBDUK step. This is the file that reads will be selected from based on their k-mer content. 
--'CIVAL' (KMC parameter) specifies the minumum number of times a k-mer must be found in each dataset to be included in the k-mer database. Seperate files will be generated for each CIVAL specified.
--'KMER_LENS' (KMC parameter) specifies the k-mer length for k-mers generated from 'samples' datasets. Seperate files are generated for each k-mer length specified. 
--'NUM_MADS:' The number of median absolute deviations away from the median a given k-mer's enrichment value much be to be included in the rest of the analysis. Separate files are generated for each NUM_MAD specified.
--'pairing_config:' specifies the two datasets that will be compared to generate enrichment values. Pairing in brackets needs to be formatted with {basename}_{readnum} (from samples parameter).
--'MEM' indicates the max memory allocated to the entire pipeline.

-run command:
--'snakemake -s {snakefile} --configfile {config_file} {options}'
--refer to snakemake documentation for options and additional parameters



