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

if [ ! -s ${outdir}/reads.fasta.fai ]
then
    samtools view ${bamfile} | awk '{ print ">"$1"\n"$10}' > ${outdir}/reads.fasta
    samtools faidx ${outdir}/reads.fasta
fi

mkdir -p ${outdir}/assembly_with_backbone

if [ ! -s ${outdir}/assembly_with_backbone/canu.contigs.fasta ]
then
    canu \
        -d ${outdir}/assembly_with_backbone \
        -p canu \
	corOutCoverage=1000 \
        minThreads=${threads} \
        useGrid=false \
        genomeSize=${genomesize} \
        gnuplotTested=true \
        -pacbio-raw ${outdir}/reads.fasta
fi

if [ ! -s ${outdir}/assembly_with_backbone/canu.contigs.quivered.fastq ]
then
    blasr \
        ${bam_file} \
        ${outdir}/assembly_with_backbone/canu.contigs.fasta \
        --bestn 1 \
        --bam \
        --nproc ${threads} \
        --out ${outdir}/assembly_with_backbone/reads_to_canu.bam
    samtools sort -@ ${threads} ${outdir}/assembly_with_backbone/reads_to_canu.bam \
        -o ${outdir}/assembly_with_backbone/reads_to_canu.sorted.bam
    pbindex ${outdir}/assembly_with_backbone/reads_to_canu.sorted.bam
    samtools faidx ${outdir}/assembly_with_backbone/canu.contigs.fasta
    quiver \
        --referenceFilename ${outdir}/assembly_with_backbone/canu.contigs.fasta \
        -j ${threads} \
        -o ${outdir}/assembly_with_backbone/canu.contigs.quivered.fasta \
        -o ${outdir}/assembly_with_backbone/canu.contigs.quivered.fastq \
        ${outdir}/assembly_with_backbone/reads_to_canu.sorted.bam
fi

if [ ! -s ${outdir}/assembly_with_backbone/backbone.bed ]
then
    sh ${ig_github}/bash/mapping.sh \
        ${backbone} \
	${outdir}/assembly_with_backbone/canu.contigs.quivered.fasta \
        1 \
	${outdir}/assembly_with_backbone/canu.contigs.quivered_to_backbone \
        10
    python ${ig_github}/python/mapping_coords.py \
        ${outdir}/assembly_with_backbone/canu.contigs.quivered_to_backbone.sorted.bam \
        > ${outdir}/assembly_with_backbone/backbone.bed 
fi

mkdir -p ${outdir}/remove_backbone/
wc -l ${outdir}/assembly_with_backbone/backbone.bed | awk '$1 == 2' | awk '{ print $2}' | while read f
do
    ref=`cat ${f} | head -1 | cut -f1`
    start=`cat ${f} | head -1 | cut -f3`
    end=`cat ${f} | tail -1 | cut -f2`
    echo ">${name}" > ${outdir}/remove_backbone/no_backbone.fasta
    samtools faidx ${outdir}/assembly_with_backbone/canu.contigs.quivered.fasta \
        ${ref}:${start}-${end} | grep -v ">" | tr -d "\n" >> \
        ${outdir}/remove_backbone/no_backbone.fasta
    echo "" >> ${outdir}/remove_backbone/no_backbone.fasta
done

wc -l ${outdir}/assembly_with_backbone/backbone.bed | awk '$1 == 1' | awk '{ print $2}' | while read f
do
    blastn \
        -query ${outdir}/assembly_with_backbone/canu.contigs.quivered.fasta \
        -subject ${outdir}/assembly_with_backbone/canu.contigs.quivered.fasta \
        -outfmt 6 > ${outdir}/remove_backbone/blast.out
    num=`cat ${outdir}/remove_backbone/blast.out | head -2 | tail -1 | cut -f9`
    if [ ${num} -lt 10 ]
    then
        ref=`cat ${outdir}/assembly_with_backbone/backbone.bed | awk '{ print $1 }'`
        start_1=`cat ${outdir}/assembly_with_backbone/backbone.bed | awk '{ print $3 }'`
        end_1=`head -2 ${outdir}/remove_backbone/blast.out | tail -1 | awk '{ print $7 }'`
        start_2=1
        end_2=`cat ${outdir}/assembly_with_backbone/backbone.bed | awk '{ print $2 }'`
        echo ">${name}" > ${outdir}/remove_backbone/no_backbone.fasta
        samtools faidx ${outdir}/assembly_with_backbone/canu.contigs.quivered.fasta \
            ${ref}:${start_1}-${end_1} | grep -v ">" | tr -d "\n" >> \
            ${outdir}/remove_backbone/no_backbone.fasta
        samtools faidx ${outdir}/assembly_with_backbone/canu.contigs.quivered.fasta \
            ${ref}:${start_2}-${end_2} | grep -v ">" | tr -d "\n" >> \
            ${outdir}/remove_backbone/no_backbone.fasta
        echo "" >> ${outdir}/remove_backbone/no_backbone.fasta
    fi
done

if [ ! -s ${outdir}/remove_backbone/no_backbone.quivered.fastq ]
then
    blasr \
        ${bam_file} \
        ${outdir}/remove_backbone/no_backbone.fasta \
        --bestn 1 \
        --bam \
        --nproc ${threads} \
        --out ${outdir}/remove_backbone/reads_to_canu.bam
    samtools sort -@ ${threads} ${outdir}/remove_backbone/reads_to_canu.bam \
        -o ${outdir}/remove_backbone/reads_to_canu.sorted.bam
    pbindex ${outdir}/remove_backbone/reads_to_canu.sorted.bam
    samtools faidx ${outdir}/remove_backbone/no_backbone.fasta
    quiver \
        --referenceFilename ${outdir}/remove_backbone/no_backbone.fasta \
        -j ${threads} \
        -o ${outdir}/remove_backbone/no_backbone.quivered.fasta \
        -o ${outdir}/remove_backbone/no_backbone.quivered.fastq \
        ${outdir}/remove_backbone/reads_to_canu.sorted.bam
fi
