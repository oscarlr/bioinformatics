#!/bin/bash

scratch=$1
work=$2
bamfile=$3
bedfile=$4
ref=$5

threads=1
memory=8

SCRIPT=$(readlink -f "$0")
SCRIPTPATH=$(dirname $(dirname "$SCRIPT"))

blasr \
    ${bamfile} \
    ${ref} \
    --bestn 1 \
    --bam \
    --nproc ${threads} \
    --out ${scratch}/reads_to_ref.bam

samtools \
    sort -@ ${threads} \
    ${scratch}/${sample}/reads_to_ref.bam \
    -o ${scratch}/${sample}/reads_to_ref.sorted.bam

samtools index ${scratch}/${sample}/reads_to_ref.sorted.bam

rm -f ${scratch}/${sample}/reads_to_ref.bam

cat ${bedfile} | while read chrom start end sample
do
    mkdir -p ${scratch}/${sample} 
    
    python ${SCRIPTPATH}/python/extract_reads_from_region.py \
	${scratch}/${sample}/reads_to_ref.sorted.bam \
        ${chrom} \
        ${start} \
        ${end} > ${scratch}/${sample}/reads.fasta
    
    samtools faidx ${scratch}/${sample}/reads.fasta
    
    # assemble fosmids

    canu \
        -d ${scratch}/${sample}/canu \
        -p canu \
        minThreads=${threads} \
        minMemory=${memory} \
        useGrid=false \
        genomeSize=50000 \
        gnuplotTested=true \
	-pacbio-raw ${scratch}/${sample}/reads.fasta
    
    if [ ! -s ${scratch}/${sample}/canu/canu.contigs.fasta ]
    then
	continue
    fi
    
    samtools faidx ${scratch}/${sample}/canu.contigs.fasta
 
    cut -f1 ${scratch}/${sample}/reads.fasta.fai | sort | uniq > ${scratch}/${sample}/reads.names
    
    python ${SCRIPTPATH}/python/extract_reads.py \
	-b ${bamfile} \
	-n ${scratch}/${sample}/reads.names \
	-o ${scratch}/${sample}/reads.bam
    
    blasr \
        ${scratch}/${sample}/reads.bam \
	${scratch}/${sample}/canu/canu.contigs.fasta \
        --bestn 1 \
        --bam \
        --nproc ${threads} \
        --out ${scratch}/${sample}/canu/reads_to_canu_contigs.bam
    
    samtools \
	sort -@ ${threads} \
	${scratch}/${sample}/canu/reads_to_canu_contigs.bam \
	-o ${scratch}/${sample}/canu/reads_to_canu_contigs.sorted.bam
    
    samtools index ${scratch}/${sample}/canu/reads_to_canu_contigs.sorted.bam
    
    rm -f ${scratch}/${sample}/canu/reads_to_canu_contigs.bam
    
    pbindex ${scratch}/${sample}/canu/reads_to_canu_contigs.sorted.bam

    quiver \
        --referenceFilename ${scratch}/${sample}/canu/canu.contigs.fasta \
        -j ${threads} \
        -o ${scratch}/${sample}/canu/canu.contigs.quiver.fasta \
        -o ${scratch}/${sample}/canu/canu.contigs.quiver.fastq \
	${scratch}/${sample}/canu/reads_to_canu_contigs.sorted.bam    
done
