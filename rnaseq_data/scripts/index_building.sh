#!/usr/bin/bash

##This script is for bulding the index of the reference genome.

module load GCC/5.4.0-2.26
module load Python/3.7.4-GCCcore-8.3.0

FOLDERGENOME=$1
ANNOTATIONDIR=$2
GENOME=$3

cd $FOLDERGENOME

extract_splice_sites.py $ANNOTATIONDIR  > splices_sites.ss
extract_exons.py $ANNOTATIONDIR  > exons.exon
hisat2-build --ss splices_sites.ss --exon exons.exon $GENOME index


