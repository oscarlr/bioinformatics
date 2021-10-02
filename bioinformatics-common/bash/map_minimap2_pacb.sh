#!/bin/bash
set -e -x

module load samtools

fastq=${1}
prefix=${2}
reffn=${3}
sample=${4}
threads=${5}

#reffn=/sc/arion/projects/sharpa01a/rodrio10/databases/references/hg38/hg38_chr_igh/minimap2_ont_index/ref.mmi

#mkdir -p ${outdir}/tmp


/sc/arion/projects/sharpa01a/rodrio10/tools/minimap2/minimap2 \
    -t ${threads} -L -a -x map-pb \
    -R "@RG\tID:${sample}\tSM:${sample}" \
    ${reffn} \
    ${fastq}  > ${prefix}.sam

samtools view -Sbh ${prefix}.sam > ${prefix}.bam

samtools sort -@ ${threads} \
	 ${prefix}.bam \
	 -o ${prefix}.sorted.bam

samtools index ${prefix}.sorted.bam

#rm -fr ${outdir}/tmp
#rm -f ${fastq_gz}
