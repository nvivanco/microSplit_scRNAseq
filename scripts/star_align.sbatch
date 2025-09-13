#!/bin/bash
#SBATCH --job-name=star_align
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --partition=normal
#SBATCH --mem=24Gb
#SBATCH --time=24:00:00

#Usage: sbatch star_align.sbatch \
#        $PROCESSED_DATA/ref_index/STAR_bsub_genome_index \
#        $RAW_DATA/sra/SRP266243/SRR11940660/SRR11940660_1.fastq \
#        $RAW_DATA/sra/SRP266243/SRR11940660/SRR11940660_2.fastq \
#        $RAW_DATA/barcodes/SRP266243/barcodes_r2_r3_solo.txt \
#        $RAW_DATA/barcodes/SRP266243/barcodes_r1_solo.txt \
#        $PROCESSED_DATA/SRP266243/SRR11940660/STAR_output


GENOME_INDEX=$1 # $PROCESSED_DATA/ref_index/STAR_bsub_genome_index
CDNA_READS=$2 # $RAW_DATA/sra/SRP266243/SRR11940660/SRR11940660_1.fastq
BC_UMI_READS=$3 # $RAW_DATA/sra/SRP266243/SRR11940660/SRR11940660_2.fastq
BC_R2_R3=$4 # $RAW_DATA/barcodes/SRP266243/barcodes_r2_r3_solo.txt
BC_R1=$5 # $RAW_DATA/barcodes/SRP266243/barcodes_r1_solo.txt
OUTPUT_DIR=$6 # $PROCESSED_DATA/SRP266243/SRR11940660/STAR_output

$RUN STAR --genomeDir $GENOME_INDEX \
          --runThreadN 12 \
          --readFilesIn $CDNA_READS $BC_UMI_READS  \
          --soloType CB_UMI_Complex \
          --soloCBposition 0_10_0_17 0_48_0_55 0_78_0_85 \
          --soloUMIposition 0_0_0_9 \
          --soloCBwhitelist $BC_R2_R3 $BC_R2_R3 $BC_R1 \
          --soloCBmatchWLtype 1MM \
          --soloUMIdedup 1MM_All \
          --soloFeatures Gene GeneFull \
          --soloMultiMappers Uniform \
          --outFilterScoreMinOverLread 0 \
          --outFilterMatchNminOverLread 0 \
          --outFilterMatchNmin 50 \
          --alignSJDBoverhangMin 1000 \
          --alignSJoverhangMin 1000 \
          --outSAMtype BAM Unsorted \
          --outFileNamePrefix $OUTPUT_DIR

#          --alignSJDBoverhangMin 1000 \ essentially disable
#          --alignSJoverhangMin 1000 \ essentially disable
