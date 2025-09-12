#!/bin/bash
#SBATCH --job-name=fastqc
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --partition=normal
#SBATCH --time=24:00:00

# Usage: sbatch fastqc.sbatch <input_folder> <output_folder>
# Example: sbatch fastqc.sbatch $RAW_DATA/sra/SRP266243/SRR11940660 $PROCESSED_DATA/SRP266243/SRR11940660/fastqc_outs

# Check for input directory
if [ $# -ne 2 ]; then
    echo "Usage: $0 <input_folder> <output_folder>"
    exit 1
fi

INPUT_DIR="$1"
OUTPUT_DIR="$2"
mkdir -p "$OUTPUT_DIR"

# Find all fastq files (compressed or not)
find "$INPUT_DIR" -type f \( -iname "*.fastq" -o -iname "*.fastq.gz" \) | while read FILE; do
    echo "Running FastQC on $FILE"
    $RUN fastqc --extract --outdir "$OUTPUT_DIR" "$FILE"
done
