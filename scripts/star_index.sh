#!/bin/bash

# Usage:  ./star_index.sh <input_fasta> <input_gtf> <output_folder>
# Example: /star_index.sh $ECOLI_REF_GENOME $ECOLI_GTF $PROCESSED_DATA/ref_index/STAR_ecoli_genome_index

# Check for input directory
if [ $# -ne 3 ]; then
    echo "Usage: $0 <input_fasta> <input_gtf> <output_folder>"
    exit 1
fi

REF_GENOME="$1"
GTF="$2"
OUTPUT_DIR="$3"

mkdir -p "$OUTPUT_DIR"

$RUN STAR --runMode genomeGenerate \
     --runThreadN 12 \
     --genomeDir $OUTPUT_DIR \
     –-genomeSAindexNbases 10 \
     --genomeFastaFiles $REF_GENOME \
     --sjdbGTFfile $GTF

### For small genomes, the parameter --genomeSAindexNbases must to be scaled down, 
### with a typical value of min(14, log2(GenomeLength)/2 - 1). 
### From https://physiology.med.cornell.edu/faculty/skrabanek/lab/angsd/lecture_notes/STARmanual.pdf