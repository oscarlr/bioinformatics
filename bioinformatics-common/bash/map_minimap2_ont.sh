#!/bin/bash
set -e -x

module load samtools

fastq=${1}
outdir=${2}
reffn=${3}
sample=${4}
threads=${5}

#reffn=/sc/arion/projects/sharpa01a/rodrio10/databases/references/hg38/hg38_chr_igh/minimap2_ont_index/ref.mmi



mkdir -p ${outdir}/tmp


/sc/arion/projects/sharpa01a/rodrio10/tools/minimap2/minimap2 \
    -t ${index} -L -a -x map-ont \
    -R "@RG\tID:${sample}\tSM:${sample}" \
    --MD \
    ${reffn} \
    ${fastq} \
    | samtools sort -@ ${index} \
	       -T ${outdir}/tmp \
	       -o ${outdir}/aln.sorted.bam

samtools index 	${outdir}/aln.sorted.bam

rm -fr ${outdir}/tmp
#rm -f ${fastq_gz}
