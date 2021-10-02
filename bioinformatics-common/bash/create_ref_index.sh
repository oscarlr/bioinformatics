#!/bin/bash
set -e -x

reference_fastafn=$1 # Full path only
outdir=$2 # Full path only

samtools faidx ${reference_fastafn}

for folder in minimap2_pb minimap2_ont pbmm2 blasr
do
    mkdir -p ${outdir}/${folder}
    if [ ! -s ${outdir}/${folder}/reference.fasta.fai ]
    then
	ln -s ${reference_fastafn} ${outdir}/${folder}/reference.fasta
	samtools faidx ${outdir}/${folder}/reference.fasta
    fi
done

minimap2 -x map-ont -d ${outdir}/minimap2_ont/reference.mmi ${outdir}/minimap2_ont/reference.fasta
minimap2 -x map-pb -d ${outdir}/minimap2_pb/reference.mmi ${outdir}/minimap2_pb/reference.fasta

sawriter ${outdir}/blasr/reference.fasta

pbmm2 index ${outdir}/pbmm2/reference.fasta ${outdir}/pbmm2/reference.mmi
