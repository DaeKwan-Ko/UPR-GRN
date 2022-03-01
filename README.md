# UPR-GRN

#The scripts provided in this repository allows to analyze raw (fastq.qz) data for RNA-seq and ChIP-seq presented in our 2022 Communications Biology paper: 

"Advanced genomics identifies growth effectors for proteotoxic ER stress recovery in Arabidopsis thaliana" (https://www.nature.com/articles/s42003-021-02964-8) 

"rna_seq/" folder contains scripts for Quality Control, Mapping, Get Read Counts and calculate differential gene expression with DESeq2.
We used the following versions of the programs:
FastQC (version 0.11.5). 
Cutadapt (version 1.8.1)
Bowtie (version 2.2.4)
TopHat (version 2.0.14)
Cufflinks (version 1.3.0)
HTSeq (version 0.6.1p1)
DESeq2 (version 1.16.1)
R (version 3.4.0)
Integrative Genome Browser (version 2.5.0)
Cytoscape software (version 3.6.1)

"chip_seq/" folder contains scripts for  Quality Control, Mapping, Peak calling and IDR
We used the following versions of the programs:
FastQC (version 0.11.5)
Cutadapt (version 1.8.1)
Bowtie (version 1.1.2)
Samtools (version 1.8)
MACS2 (version 2.1.2)
R (version 3.4.0)
ChIPseeker (version 3.12)
