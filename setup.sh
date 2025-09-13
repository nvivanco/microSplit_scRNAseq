#!/bin/bash
# default paths
export RAW_DATA="$OAK/scRNAseq/raw_data/"
export PROCESSED_DATA="$OAK/scRNAseq/processed_data/"
export ECOLI_REF_GENOME="$RAW_DATA/references/ecoli/GCF_000005845.2_ASM584v2_genomic.fna"
export ECOLI_GTF="$RAW_DATA/references/ecoli/GCF_000005845.2_ASM584v2_genomic.gtf"
export ECOLI_REF_TRNX="$RAW_DATA/references/ecoli/GCF_000005845.2_ASM584v2_cds_from_genomic.fna"

export BSUB_REF_GENOME="$RAW_DATA/references/bsubtilis/GCF_000009045.1_ASM904v1_genomic.fna"
export BSUB_GTF="$RAW_DATA/references/bsubtilis/GCF_000009045.1_ASM904v1_genomic.gtf"
export BSUB_REF_TRNX="$RAW_DATA/references/bsubtilis/GCF_000005845.GCF_000009045.1_ASM904v1_cds_from_genomic.fna"

# run with singularity
export RUN="singularity run $OAK/scRNAseq/envs/scrna_latest.sif"