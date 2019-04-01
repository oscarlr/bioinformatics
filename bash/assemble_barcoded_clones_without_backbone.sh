#!/bin/bash
set -x -e

outdir=$1
bamfile=$2
genomesize=$3
backbone=$4
name=$5

threads=20

ig_github=/sc/orga/work/rodrio10/software/in_github/ig

mkdir -p ${outdir}

if [ ! -s ${outdir}/reads_to_backbone.sorted.bam.bai ]
then
    blasr \
	${bamfile} \
	${backbone} \
	--unaligned ${outdir}/reads_to_backbone.unaligned.fasta \
	--bestn 1 \
	--bam \
	--nproc ${threads} \
	--out ${outdir}/reads_to_backbone.bam
    samtools sort -@ ${threads} ${outdir}/reads_to_backbone.bam \
	-o ${outdir}/reads_to_backbone.sorted.bam
    samtools index ${outdir}/reads_to_backbone.sorted.bam
    rm -f ${threads} ${outdir}/reads_to_backbone.bam
fi

if [ ! -s ${outdir}/reads.fasta ]
then
    cat ${outdir}/reads_to_backbone.unaligned.fasta > ${outdir}/reads.fasta
    python ${ig_github}/python/remove_backbone.py \
	${outdir}/reads_to_backbone.sorted.bam >> ${outdir}/reads.fasta
fi

if [ ! -s ${outdir}/assembly_without_backbone/canu.contigs.fasta ]
then
    canu \
        -d ${outdir}/assembly_without_backbone \
        -p canu \
	corOutCoverage=200 \
        minThreads=${threads} \
        useGrid=false \
        genomeSize=${genomesize} \
        gnuplotTested=true \
        -pacbio-raw ${outdir}/reads.fasta
fi

if [ ! -s ${outdir}/assembly_without_backbone/canu.contigs.quivered.fastq ]
then
    blasr \
        ${bamfile} \
        ${outdir}/assembly_without_backbone/canu.contigs.fasta \
        --bestn 1 \
        --bam \
        --nproc ${threads} \
        --out ${outdir}/assembly_without_backbone/reads_to_canu.bam
    samtools sort -@ ${threads} ${outdir}/assembly_without_backbone/reads_to_canu.bam \
        -o ${outdir}/assembly_without_backbone/reads_to_canu.sorted.bam
    pbindex ${outdir}/assembly_without_backbone/reads_to_canu.sorted.bam
    samtools faidx ${outdir}/assembly_without_backbone/canu.contigs.fasta
    quiver \
        --referenceFilename ${outdir}/assembly_without_backbone/canu.contigs.fasta \
        -j ${threads} \
        -o ${outdir}/assembly_without_backbone/canu.contigs.quivered.fasta \
        -o ${outdir}/assembly_without_backbone/canu.contigs.quivered.fastq \
        ${outdir}/assembly_without_backbone/reads_to_canu.sorted.bam
    rm -f ${outdir}/assembly_without_backbone/reads_to_canu.bam
fi
