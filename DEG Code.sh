#!/usr/bin/bash

##########################################################
# Transcriptome Analysis Script for Diffenrial gene regulation Research
# Author: Vikas
# Date: 22/11/2023
# Description: This script performs DGE analyses including FastQC, HISAT2 alignment, stringtie assembly and abundance  
# for comparison between male and female samples from two geographcal distict population.
##########################################################

# 1. Bash script for running FastQC on multiple samples, FastQC generates HTML reports for easy interpretation and visual inspection.

SAMPLES="ERR188044 ERR188104 ERR188234 ERR188245 ERR188257 ERR188273 ERR188337 ERR188383 ERR188401 ERR188428 ERR188454 ERR204916"

for SAMPLE in $SAMPLES; do
    fastqc chrX_data/samples/${SAMPLE}_chrX_1.fastq.gz chrX_data/samples/${SAMPLE}_chrX_2.fastq.gz -o Fastqc/
done

# 2. Bash script for HISAT2; aligning all .fastq.gz files to the indexed reference genome to generate .sam files

for SAMPLE in $SAMPLES; do
    hisat2 -p 8 --dta -x chrX_data/indexes/chrX_tran -1 chrX_data/samples/${SAMPLE}_chrX_1.fastq.gz -2 chrX_data/samples/${SAMPLE}_chrX_2.fastq.gz -S map/${SAMPLE}_chrX.sam
done

# 3. Bash script for SAMTools; converting .sam files to .bam files

for SAMPLE in $SAMPLES; do
    samtools sort -@ 8 -o map/${SAMPLE}_chrX.bam map/${SAMPLE}_chrX.sam
done

# 4. Bash script for StringTie; creating assembly per sample

for SAMPLE in $SAMPLES; do
    stringtie map/${SAMPLE}_chrX.bam -l ERR188044 -p 1 -G chrX_data/genes/chrX.gtf -o assembly/${SAMPLE}_chrX.gtf
done

# 5. Bash script for getting abundance using StringTie

for SAMPLE in $SAMPLES; do
    stringtie -e -B -p 8 -G stringtie_merged.gtf -o ballgown/${SAMPLE}/${SAMPLE}_chrX.gtf map/${SAMPLE}_chrX.bam
done
