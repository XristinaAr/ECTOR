#!/usr/bin/bash

module load GCC/5.4.0-2.26
module load Python/3.7.4-GCCcore-8.3.0

FOLDER=$1
SAMPLE=$2

## Sample quality control and read mapping to reference genome

cd $FOLDER

/home/arvanitidou/.local/bin/FastQC/fastqc $SAMPLE.fastq.gz

hisat2 --dta -x /home/arvanitidou/data/temp_ostta_rnaseq/genome/index -U $SAMPLE.fastq.gz -S $SAMPLE.sam

## Generting sorted bam file

/home/arvanitidou/opt/samtools sort -o $SAMPLE.bam $SAMPLE.sam
/home/arvanitidou/opt/samtools index $SAMPLE.bam

bamCoverage -bs 5 --normalizeUsing CPM --bam $SAMPLE.bam -o $SAMPLE.bw

## Transcript assembly
/home/arvanitidou/.local/bin/stringtie-2.2.1.Linux_x86_64/stringtie -G /home/arvanitidou/data/temp_ostta_rnaseq/annotation/ostreococcus_tauri.gtf -o $SAMPLE.gtf -l $SAMPLE $SAMPLE.bam

## Gene Expression Quantification
/home/arvanitidou/.local/bin/stringtie-2.2.1.Linux_x86_64/stringtie -e -B -G /home/arvanitidou/data/temp_ostta_rnaseq/annotation/ostreococcus_tauri.gtf -o $SAMPLE.gtf $SAMPLE.bam


