#!/bin/bash
#module purge
#module load python py_packages samtools
set -e -x

reads=$1
ref=$2
threads=$3
output=$4
n=$5


if [ ! -s ${output}.sorted.bam.bai ]
then
    blasr \
    	--bestn ${n} \
    	--clipping soft \
    	--sam \
    	--nproc ${threads} \
    	--noSplitSubreads \
    	--unaligned ${output}.unaligned.fasta \
    	${reads} \
    	${ref} \
    	--out ${output}.sam
    
    samtools view -F 4 -Sbh ${output}.sam > ${output}.bam 
    samtools sort ${output}.bam -o  ${output}.sorted.bam
    samtools index ${output}.sorted.bam
    
    rm -f ${output}.sam
    rm -f ${output}.bam
fi