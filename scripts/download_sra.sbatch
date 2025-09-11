#!/bin/bash
#SBATCH --job-name=download_sra
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --partition=normal
#SBATCH --time=2:00:00

#Example usage: ./download_sra.sh SRP SRP099835 $SCRATCH/data/sra
#$1: file_type (SRA or SRP)
#$2: accession number
#$3: outputdir

ml biology load sra-tools
# Check if an argument is provided
if [ $# -eq 0 ]; then
    echo "Error: No argument provided."
    echo "Usage: $0 <argument>"
    exit 1
fi
# Store the argument in a variable
file_type="$1"
accession="$2"
output_dir="$3"
# Perform conditional based on the value of the argument
if [ "$file_type" = "SRP" ]; then

    mkdir "$output_dir/$accession"
    echo "Downloading SRP accession: $2"
    prefetch $accession -O $output_dir/$accession
    for file in $output_dir/$accession/*; do
        echo "$file"
        fasterq-dump -O $file $file
    done
elif [ "$file_type" = "SRA" ]; then
        fasterq-dump -O $output_dir $accession
else
    echo "Invalid argument: $argument"
    echo "Usage: $0 <option1|option2>"
    exit 1
fi
