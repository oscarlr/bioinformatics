# Fosmids
## Requirements
1. Python version 2
2. BLAST/2.7.1+
3. WhatsHap: Look at https://github.com/oscarlr/IGenotyper <br>
4. BLASR

## Code to merge the fosmids
```
# fosmids.fasta has all the fosmids sequence to be merged
module load blast/2.7.1+
blastn \
  -query fosmids.fasta \
  -subject fosmids.fasta \
  -outfmt "6 length pident nident mismatch gapopen gaps qseqid qstart qend qlen sseqid sstart send slen sstrand" > \
  blast.txt
 
# blast.txt is the output from the blast output
# blast_edited.txt is the blast
# fosmid_groups.txt are a txt file with fosmids that belong together
# fosmids_to_ignore.txt are a txt file with fosmids that should be ignored
# fosmids_to_merge.txt are a txt file with fosmids that are being merged
# fosmids.fasta is the input fasta
# merged_fosmids.fasta is the output fasta
# 5000 the minimum number of bases to overlap
# 1 is the maximum allowed errors
python python/merge_fosmids.py 
  blast.txt \
  blast_edited.txt \
  fosmid_groups.txt \
  fosmids_to_ignore.txt \
  fosmids_to_merge.txt \
  fosmids.fasta \
  merged_fosmids.fasta \
  5000 \
  1
 
# fosmids_to_ref.sorted.bam is the bam file with the fosmids aligned
# fosmids_to_ref_grouped.sorted.bam is the bam file with the fosmids that are to be merged
# fosmid_groups.txt is the output from above
python python/add_read_group.py \
  fosmids_to_ref.sorted.bam \
  fosmids_to_ref_grouped.sorted.bam \
  fosmid_groups.txt
samtools index fosmids_to_ref_grouped.sorted.bam

```
## Aligning merged fosmids to reference
```
# threads is the number of threads to use
# this will create a file "merged_fosmids_to_ref.sorted.bam" containing the merged fosmids aligned to the reference
sh bash/mapping_indel_penalize.sh \
   merged_fosmids.fasta \
   ref.fasta \
   ${threads} \
   merged_fosmids_to_ref \
   1
```
## Extracting phasing SNVs from 1KG
```
tmp=
mkdir -p ${tmp}/get_1kg_snvs
wget -O ${tmp}/get_1kg_snvs/ALL.chr14_GRCh38.genotypes.20170504.vcf.gz \
  http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/GRCh38_positions/ALL.chr14_GRCh38.genotypes.20170504.vcf.gz
vcftools \
  --gzvcf ${tmp}/get_1kg_snvs/ALL.chr14_GRCh38.genotypes.20170504.vcf.gz \
  --chr 14 \
  --from-bp 105917017 \
  --to-bp 107044718 \
  --indv NA18507 \
  --recode --recode-INFO-all --stdout \
  | grep -v "0|0" | grep -v "1|1" | grep -v "./." | sed 's/^14/chr14/g' > ${tmp}/get_1kg_snvs/NA18507.vcf
```

## Rephasing SNVs with fosmids
```
tmp=
mkdir -p ${tmp}/rephase_1kg_snvs
ref=hg38.fa
fosmids=fosmids_to_hg38.sorted.bam
vcf=${tmp}/get_1kg_snvs/NA18507.vcf
whatshap phase \
    --sample NA18507 \
    --reference ${ref} \
    --ignore-read-groups \
    --distrust-genotypes \
    -o ${tmp}/snvs_in_hg38/snps_phased.vcf \
    ${tmp}/snvs_in_hg38/snps.vcf \
    ${vcf} \
    ${fosmids}
```

## Phasing fosmids with rephased SNVs
```
tmp=
mkdir -p ${tmp}/rephase_1kg_snvs_reads
fosmids=fosmids_to_hg38.sorted.bam
python python/eval_phase_het_snps.py \
   ${tmp}/snvs_in_hg38/snps_phased.vcf \
   ${fosmids} \
   ${tmp}/rephase_1kg_snvs_reads/fosmids_to_hg38_phased.sorted.bam \
   NA18507 > ${tmp}/rephase_1kg_snvs_reads/fosmids_to_hg38_phased_scores.bed
samtools index ${tmp}/rephase_1kg_snvs_reads/fosmids_to_hg38_phased.sorted.bam

```
