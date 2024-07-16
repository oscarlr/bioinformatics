#!/bin/bash
set -e -x

fastq=${1}
prefix=${2}
reffn=${3}
sample=${4}
threads=${5}

#reffn=/sc/arion/projects/sharpa01a/rodrio10/databases/references/hg38/hg38_chr_igh/minimap2_ont_index/ref.mmi

#mkdir -p ${outdir}/tmp

# map-pac


base=`echo ${prefix} | rev | cut -f2- -d/ | rev`

tmpdir="${base}/$RANDOM"

#    -O 2,15 \

#    -w 5 \
#    -O 2,15 \

/home/o0rodr03/downloads/minimap2/minimap2 \
    -Y \
    --secondary=no \
    --split-prefix ${tmpdir} \
    -x asm20 \
    -t ${threads} -L -a \
    -R "@RG\tID:${sample}\tSM:${sample}" \
    ${reffn} \
    ${fastq}  > ${prefix}.sam

samtools view -Sbh ${prefix}.sam > ${prefix}.bam

samtools sort -@ ${threads} \
	 ${prefix}.bam \
	 -o ${prefix}.sorted.bam

samtools index ${prefix}.sorted.bam

rm -f ${prefix}.sam
rm -f ${prefix}.bam
rm -fr ${tmpdir}

#rm -fr ${outdir}/tmp
#rm -f ${fastq_gz}
