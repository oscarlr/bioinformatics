#!/bin/bash
set -e -x

scratch=$1
bamfile=$2
chrom=$3
start=$4
end=$5
sample=$6
alignment=$7

threads=1
memory=8

SCRIPT=$(readlink -f "$0" )
SCRIPTPATH=$(dirname $(dirname "$SCRIPT"))

mkdir -p ${scratch}/${sample} 

if [ ! -s ${scratch}/${sample}/canu/canu.contigs.fasta  ]
then
    python ${SCRIPTPATH}/python/extract_reads_from_region.py \
	${alignment} \
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
	useGrid=0 \
	stopOnLowCoverage=0 \
	genomeSize=50000 \
	-pacbio-raw ${scratch}/${sample}/reads.fasta
fi

if [ ! -s ${scratch}/${sample}/canu/canu.contigs.fasta ]
then
    continue
fi

if [ ! -s ${scratch}/${sample}/canu/canu.contigs.quiver.fastq ]
then
    samtools faidx ${scratch}/${sample}/canu/canu.contigs.fasta
    
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

    samtools faidx ${scratch}/${sample}/canu/canu.contigs.quiver.fasta
fi

if [ ! -s ${scratch}/${sample}/${sample}.fasta ]
then
    
    sh ${SCRIPTPATH}/bash/mapping.sh \
	/sc/orga/work/rodrio10/bashia02c_projects/igh/humancomplexityigh/permanent/fosmid_backbone_pcc2_site.fa \
	${scratch}/${sample}/canu/canu.contigs.quiver.fasta \
	1 \
	${scratch}/${sample}/canu/backbone_to_canu.contigs.quiver \
	10
    
    python ${SCRIPTPATH}/python/mapping_coords.py \
	${scratch}/${sample}/canu/backbone_to_canu.contigs.quiver.sorted.bam \
	> ${scratch}/${sample}/canu/backbone.bed
    
    wc -l ${scratch}/${sample}/canu/backbone.bed  | awk '$1 == 2' | awk '{ print $2}' | while read f
    do
	ref=`cat ${f} | head -1 | cut -f1`
	start=`cat ${f} | head -1 | cut -f3`
	end=`cat ${f} | tail -1 | cut -f2`
	echo ">${sample}" > ${scratch}/${sample}/${sample}.fasta
	samtools faidx ${scratch}/${sample}/canu/canu.contigs.quiver.fasta \
            ${ref}:${start}-${end} | grep -v ">" | tr -d "\n" >> \
            ${scratch}/${sample}/${sample}.fasta
	echo "" >> ${scratch}/${sample}/${sample}.fasta
    done
    
    wc -l ${scratch}/${sample}/canu/backbone.bed  | awk '$1 == 1' | awk '{ print $2}' | while read f
    do
	blastn \
            -query ${scratch}/${sample}/canu/canu.contigs.quiver.fasta \
            -subject ${scratch}/${sample}/canu/canu.contigs.quiver.fasta \
            -outfmt 6 > ${scratch}/${sample}/blast.out
	num=`cat ${scratch}/${sample}/blast.out | head -2 | tail -1 | cut -f9`
	if [ ${num} -lt 10 ]
	then
            ref=`cat ${scratch}/${sample}/canu/backbone.bed | awk '{ print $1 }'`
            start_1=`cat ${scratch}/${sample}/canu/backbone.bed | awk '{ print $3 }'`
            end_1=`head -2 ${scratch}/${sample}/blast.out | tail -1 | awk '{ print $7 }'`
            start_2=1
            end_2=`cat ${scratch}/${sample}/canu/backbone.bed | awk '{ print $2 }'`
            echo ">${sample}" > ${scratch}/${sample}/${sample}.fasta
            samtools faidx ${scratch}/${sample}/canu/canu.contigs.quiver.fasta \
		${ref}:${start_1}-${end_1} | grep -v ">" | tr -d "\n" >> \
		${scratch}/${sample}/${sample}.fasta
            samtools faidx ${scratch}/${sample}/canu/canu.contigs.quiver.fasta \
		${ref}:${start_2}-${end_2} | grep -v ">" | tr -d "\n" >> \
		${scratch}/${sample}/${sample}.fasta
            echo "" >> ${scratch}/${sample}/${sample}.fasta
	else
	    ref=`cat ${scratch}/${sample}/canu/backbone.bed | awk '{ print $1 }'`
            backbone_start=`cat ${scratch}/${sample}/canu/backbone.bed | awk '{ print $2 }'`
            backbone_end=`cat ${scratch}/${sample}/canu/backbone.bed | awk '{ print $3 }'`
	    contig_length=`cut -f2 ${scratch}/${sample}/canu/canu.contigs.quiver.fasta.fai`
	    distance_to_end=`echo -e "${backbone_end}\t${contig_length}" | awk '{ print $2 - $1}'`
	    if [ ${distance_to_end} -lt ${backbone_start} ]
	    then
		echo ">${sample}" > ${scratch}/${sample}/${sample}.fasta
		samtools faidx ${scratch}/${sample}/canu/canu.contigs.quiver.fasta \
		    ${ref}:1-${backbone_start} | grep -v ">" | tr -d "\n" >> \
		    ${scratch}/${sample}/${sample}.fasta
	    else
		echo ">${sample}" > ${scratch}/${sample}/${sample}.fasta
		samtools faidx ${scratch}/${sample}/canu/canu.contigs.quiver.fasta \
		    ${ref}:${backbone_end}-${contig_length} | grep -v ">" | tr -d "\n" >> \
		    ${scratch}/${sample}/${sample}.fasta		
	    fi
	fi
    done
fi

if [ ${sample} == "ABC9_3_2_000043875600_C12" ]
then
    echo ">${sample}" > ${scratch}/${sample}/${sample}.fasta
    samtools faidx ${scratch}/${sample}/canu/canu.contigs.quiver.fasta "tig00000007|quiver":1845-44799 >> \
	${scratch}/${sample}/${sample}.fasta
fi

if [ ${sample} == "ABC10_2_1_000044377500_N2" ]
then
    echo ">${sample}" > ${scratch}/${sample}/${sample}.fasta
    samtools faidx ${scratch}/${sample}/canu/canu.contigs.quiver.fasta "tig00000001|quiver":16312-55089 >> \
	${scratch}/${sample}/${sample}.fasta
fi

if [ ${sample} == "1205562_ABC14_3_1_000001180022_A15" ]
then
    echo ">${sample}" > ${scratch}/${sample}/${sample}.fasta
    samtools faidx ${scratch}/${sample}/canu/canu.contigs.quiver.fasta "tig00000010|quiver":11234-52056 >> \
	${scratch}/${sample}/${sample}.fasta
fi

if [ ! -s ${scratch}/${sample}/${sample}.quiver.fastq ]
then
    blasr \
	${scratch}/${sample}/reads.bam \
	${scratch}/${sample}/${sample}.fasta \
	--bestn 1 \
	--bam \
	--nproc ${threads} \
	--out ${scratch}/${sample}/reads_to_${sample}.bam
    
    samtools \
	sort -@ ${threads} \
	${scratch}/${sample}/reads_to_${sample}.bam \
	-o ${scratch}/${sample}/reads_to_${sample}.sorted.bam
    
    samtools index ${scratch}/${sample}/reads_to_${sample}.sorted.bam
    
    rm -f ${scratch}/${sample}/reads_to_${sample}.bam
    
    pbindex ${scratch}/${sample}/reads_to_${sample}.sorted.bam
    
    samtools faidx ${scratch}/${sample}/${sample}.fasta
    
    quiver \
	--referenceFilename ${scratch}/${sample}/${sample}.fasta \
	-j ${threads} \
	-o ${scratch}/${sample}/${sample}.quiver.fasta \
	-o ${scratch}/${sample}/${sample}.quiver.fastq \
	${scratch}/${sample}/reads_to_${sample}.sorted.bam    
fi
