#!/bin/bash
set -e -x

reads=$1
ref=$2
threads=$3
output=$4
n=$5


if [ ! -s ${output}.sorted.bam.bai ]
then
    blasr \
	--sa ${ref}.sa \
    	--nproc ${threads} \
    	--unaligned ${output}.unaligned.fasta \
    	--out ${output}.sam \
    	${reads} \
    	${ref} \
	--noSplitSubreads \
	--insertion 0 --deletion 0 \
	--minMatch 15 \
	--maxMatch 30 \
	--advanceHalf \
	--fastSDP \
        --sam 
    
    samtools view -F 4 -Sbh ${output}.sam > ${output}.bam 
    samtools sort ${output}.bam -o  ${output}.sorted.bam
    samtools index ${output}.sorted.bam
    
    rm -f ${output}.sam
    rm -f ${output}.bam
fi
