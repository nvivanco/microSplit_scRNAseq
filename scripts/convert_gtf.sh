#!/bin/bash
input_gtf=$1
output_gtf="${input_gtf%.gtf}.converted.gtf"
awk 'BEGIN{OFS="\t"} $3=="CDS"{$3="exon"} {print}' $input_gtf > $output_gtf