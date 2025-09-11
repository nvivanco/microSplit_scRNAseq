#!/bin/bash

# Usage: ./download_refgenome.sh <output_directory_path>
# Example: ./download_refgenome.sh /data/my_project

# Check if an argument was provided
if [ -z "$1" ]; then
    echo "Usage: $0 <output_directory_path>"
    exit 1
fi

# Assign the first argument to a variable
output_dir="$1"

# Create output directories using the variable
mkdir -p "$output_dir/references/ecoli"
mkdir -p "$output_dir/references/bsubtilis"

echo "Directories created successfully in $output_dir"

# --- Download E. coli Reference Files ---
echo "Downloading E. coli reference files..."
wget -P "$output_dir/references/ecoli" https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz
wget -P "$output_dir/references/ecoli" https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.gtf.gz
wget -P "$output_dir/references/ecoli" https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_cds_from_genomic.fna.gz

# --- Download B. subtilis Reference Files ---
echo "Downloading B. subtilis reference files..."
wget -P "$output_dir/references/bsubtilis" https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/009/045/GCF_000009045.1_ASM904v1/GCF_000009045.1_ASM904v1_genomic.fna.gz
wget -P "$output_dir/references/bsubtilis" https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/009/045/GCF_000009045.1_ASM904v1/GCF_000009045.1_ASM904v1_genomic.gtf.gz
wget -P "$output_dir/references/bsubtilis" https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/009/045/GCF_000009045.1_ASM904v1/GCF_000009045.1_ASM904v1_cds_from_genomic.fna.gz

# --- Unzip all downloaded files ---
echo "Unzipping downloaded files..."
gunzip "$output_dir/references/ecoli"/*.gz
gunzip "$output_dir/references/bsubtilis"/*.gz

echo "Script finished."