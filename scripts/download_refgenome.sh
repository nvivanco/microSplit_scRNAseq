#!/bin/bash

# Usage: ./download_refgenome 
# Run this script in directory where the output files should be

# Create output directory
mkdir -p references/ecoli
mkdir -p references/bsubtilis

# --- Download E. coli Reference Files ---
wget -P references/ecoli https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz
wget -P references/ecoli https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.gtf.gz
wget -P references/ecoli https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_cds_from_genomic.fna.gz

# --- Download B. subtilis Reference Files ---
wget -P references/bsubtilis https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/009/045/GCF_000009045.1_ASM904v1/GCF_000009045.1_ASM904v1_genomic.fna.gz
wget -P references/bsubtilis https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/009/045/GCF_000009045.1_ASM904v1/GCF_000009045.1_ASM904v1_genomic.gtf.gz
wget -P references/bsubtilis https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/009/045/GCF_000009045.1_ASM904v1/GCF_000009045.1_ASM904v1_cds_from_genomic.fna.gz

# --- Unzip all downloaded files ---
gunzip references/ecoli/*.gz
gunzip references/bsubtilis/*.gz
