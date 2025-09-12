#!/bin/bash
# default paths
export RAW_DATA="$OAK/scRNAseq/raw_data/"
export PROCESSED_DATA="$OAK/scRNAseq/processed_data/"

# run with singularity
export RUN="singularity run $OAK/scRNAseq/envs/scrna_latest.sif"