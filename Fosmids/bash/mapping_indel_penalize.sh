#!/bin/bash
#module purge
#module load python py_packages samtools
set -e -x

reads=$1
ref=$2
threads=$3
output=$4
n=$5

#BLASR_PATH=/hpc/users/bashia02/gitrepos/mchaisson7/blasr/alignment/bin # previous
#BLASR_PATH=/sc/orga/work/rodrio10/tools/blasr/alignment/bin


    # blasr \
    # 	--nproc ${threads} \
    # 	--unaligned ${output}.unaligned.fasta \
    # 	--out ${output}.sam \
    # 	${reads} \
    # 	${ref} \
    #     --noSplitSubreads \
    #     --sam \
    #     --nproc 1 \
    #     --minMatch 5 \
    #     --maxMatch 20 \
    #     --advanceHalf \
    #     --advanceExactMatches 10 \
    #     --fastMaxInterval \
    #     --fastSDP \
    #     --aggressiveIntervalCut

if [ ! -s ${output}.sorted.bam.baix ]
then
    # --maxMatch 15 \
    blasr \
	--sa ${ref}.sa \
    	--nproc ${threads} \
    	--unaligned ${output}.unaligned.fasta \
    	--out ${output}.sam \
    	${reads} \
    	${ref} \
	--noSplitSubreads \
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
