#!/bin/bash
set -e -x

sample=$1
reffn=$2
read_alignment=$3
chrom=$4
outdir=$5

tech="--pacbio"

whatshap find_snv_candidates \
    --sample ${sample} \
    ${reffn} \
    ${read_alignment} \
    ${tech} \
    --chromosome "${chrom}" \
    -o ${outdir}/snv_candidates.vcf

whatshap genotype \
	 --sample ${sample} \
	 --ignore-read-groups \
	 --reference ${reffn} \
	 -o ${outdir}/snvs.vcf \
	 ${outdir}/snv_candidates.vcf \
	 ${read_alignment}
