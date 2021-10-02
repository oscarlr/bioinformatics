#!/bin/bash

scratch=$1

fosmid=$2
backbone=$3
fosmid_name=$4
reads_bam_format=$5

threads=1

SCRIPT=$(readlink -f "$0")
SCRIPTPATH=$(dirname $(dirname "$SCRIPT"))

blasr \
    ${backbone} \
    ${fosmid} \
    --bestn 10 \
    --bam \
    --nproc ${threads} \
    --out ${scratch}/backbone_to_fosmid.bam

samtools \
    sort -@ ${threads} \
    ${scratch}/backbone_to_fosmid.bam \
    -o ${scratch}/backbone_to_fosmid.sorted.bam

samtools index ${scratch}/backbone_to_fosmid.sorted.bam

rm -f ${scratch}/backbone_to_fosmid.bam

python ${SCRIPTPATH}/mapping_coords.py \
    ${scratch}/backbone_to_fosmid.sorted.bam > ${scratch}/backbone_to_fosmid.bed

wc -l ${scratch}/backbone_to_fosmid.bed | awk '$1 == 2' | awk '{ print $2 }' | while read f
do
    ref=`cat ${f} | head -1 | cut -f1`
    start=`cat ${f} | head -1 | cut -f3`
    end=`cat ${f} | tail -1 | cut -f2`
    echo ">${fosmid_name}" > ${scratch}/${fosmid_name}_no_backbone.fasta 
    samtools faidx ${fosmid} ${ref}:${start}-${end} | grep -v ">" | tr -d "\n" >> ${scratch}/${fosmid_name}_no_backbone.fasta 
    echo "" >> ${scratch}/${fosmid_name}_no_backbone.fasta 
done

wc -l ${scratch}/backbone_to_fosmid.bed | awk '$1 == 1' | awk '{ print $2 }' | while read f
do
    blastn \
	-query ${fosmid} \
	-subject ${fosmid} \
	-outfmt 6 > ${scratch}/${fosmid_name}.blast
    
    num=`cat ${scratch}/${fosmid_name}.blast | head -2 | tail -1 | cut -f9`
    if [ ${num} -lt 10 ]
    then
        ref=`cat ${scratch}/backbone_to_fosmid.bed | awk '{ print $1 }'`
        start_1=`cat ${scratch}/backbone_to_fosmid.bed | awk '{ print $3 }'`
        end_1=`head -2 ${scratch}/${fosmid_name}.blast | tail -1 | awk '{ print $7 }'`
        start_2=1
        end_2=`cat ${scratch}/backbone_to_fosmid.bed | awk '{ print $2 }'`
        echo ">${fosmid_name}" > ${scratch}/${fosmid_name}_no_backbone.fasta 
        samtools faidx ${fosmid} \
            ${ref}:${start_1}-${end_1} | grep -v ">" | tr -d "\n" >> \
            ${scratch}/${fosmid_name}_no_backbone.fasta 
        samtools faidx ${fosmid} \
            ${ref}:${start_2}-${end_2} | grep -v ">" | tr -d "\n" >> \
            ${scratch}/${fosmid_name}_no_backbone.fasta 
        echo "" >> ${scratch}/${fosmid_name}_no_backbone.fasta
    fi
done


if [ -s ${scratch}/${fosmid_name}_no_backbone.fasta ]
then
    samtools faidx ${scratch}/${fosmid_name}_no_backbone.fasta 

    blasr \
	${reads_bam_format} \
	${scratch}/${fosmid_name}_no_backbone.fasta \
        --bestn 1 \
        --bam \
        --nproc ${threads} \
        --out ${scratch}/reads_to_${fosmid_name}_no_backbone.bam
    
    samtools \
	sort -@ ${threads} \
	${scratch}/reads_to_${fosmid_name}_no_backbone.bam \
	-o ${scratch}/reads_to_${fosmid_name}_no_backbone.sorted.bam
    
    samtools index ${scratch}/reads_to_${fosmid_name}_no_backbone.sorted.bam
    
    rm -f ${scratch}/reads_to_${fosmid_name}_no_backbone.bam
    
    pbindex ${scratch}/reads_to_${fosmid_name}_no_backbone.sorted.bam

    quiver \
        --referenceFilename ${scratch}/${fosmid_name}_no_backbone.fasta \
        -j ${threads} \
        -o ${scratch}/${fosmid_name}_no_backbone.quiver.fasta \
        -o ${scratch}/${fosmid_name}_no_backbone.quiver.fastq \
	${scratch}/reads_to_${fosmid_name}_no_backbone.sorted.bam  
    
fi

