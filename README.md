# bioinformatics
## Mergining contigs

This code is useful to merge contigs (generally any sequences). We have used this code to merge contigs, align back to the reference and detect SVs.

The steps are:
1. Perform a self-alignment using BLAST
2. Merge the sequences using the BLAST output (merge_fosmids.py)

The sequences selected to be be merged can be evaluated using add_read_group.py


# fosmids.fasta has all the fosmids sequence to be merged
module load blast/2.7.1+
blastn \
  -query seqs.fasta \
  -subject seqs.fasta \
  -outfmt "6 length pident nident mismatch gapopen gaps qseqid qstart qend qlen sseqid sstart send slen sstrand" > \
  blast.txt
 
# blast.txt is the output from the blast output
# blast_edited.txt is the blast
# seqs_groups.txt are a txt file with fosmids that belong together
# seqs_to_ignore.txt are a txt file with fosmids that should be ignored
# seqs_to_merge.txt are a txt file with fosmids that are being merged
# seqs.fasta is the input fasta
# merged_seqs.fasta is the output fasta
# 5000 the minimum number of bases to overlap
# 1 is the maximum allowed errors
python Fosmids/python/merge_fosmids.py 
  blast.txt \
  blast_edited.txt \
  seqs_groups.txt \
  seqs_to_ignore.txt \
  seqs_to_merge.txt \
  seqs.fasta \
  merged_seqs.fasta \
  5000 \
  1
 
# seqs_to_ref.sorted.bam is the bam file with the fosmids aligned
# seqs_to_ref_grouped.sorted.bam is the bam file with the fosmids that are to be merged
# seq_groups.txt is the output from above
python Fosmids/python/add_read_group.py \
  seqs_to_ref.sorted.bam \
  seqs_to_ref_grouped.sorted.bam \
  seqs_groups.txt
samtools index seqs_to_ref_grouped.sorted.bam



### Parse IGenotyper alleles using IMGT V-Quest
```
python validate_allele_calls.py -h
usage: validate_allele_calls.py [-h] [--species SPECIES] [--locus LOCUS]
                                assembly_genes ccs_genes samples_fofn
                                airr_vquest_output

Process IGenotyper alleles

positional arguments:
  assembly_genes      assembly_genes.fasta IGenotyper output
  ccs_genes           ccs_genes.fasta IGenotyper output
  samples_fofn        A single column file with file paths of
                      assemble_genes.fasta for other samples
  airr_vquest_output  AIRR IMGT v-quest output file

optional arguments:
  -h, --help          show this help message and exit
  --species SPECIES   Species. Must be an IMGT option
  --locus LOCUS       Receptor or Locus type. Must be an IMGT option
```


### Get sequence (fasta file) from bam file for each bed entry
`python/extract_sequence_from_bam.py <bed file> <bam file> > <fasta file>`

### Get tab delimited file with number of reads per bed entry file
`python/check_mapping_in_coords.py <bed file> <bam file> > <tab-delim file>`

### Get regions (bed file) with coverage
`python/get_pos_with_bases.py <bam file> > <bed file>`

### Calculate true positive, false positeve rate, and number of false positive and true negative SNVs
#### LENGTH is the amount of bases compared. For example, LENGTH is the number of bases that overlap between two assemblies.
`python/compare_test_vcf_to_true_vcf.py <true VCF> <test VCF> <LENGTH>`

### Extract reads by read name from bam file (code found online)
`python/extract_reads.py -b <bam file> -n <list of read names> -o <bam file>`

### Extract sequences from region
`python/extract_reads_from_region.py <bam file> <chrom> <start> <end> > <fasta file>`

### Get mapping coordinates of bam entries
`python/mapping_coords.py <bam file> > <bed file>`


## bash scripts

### Create bam file with BLASR
`bash/mapping.sh <reads fasta file> <ref fasta file> <threads> <prefix> <n>`

### Assemble fosmids
`bash/assemble_fosmids.sh <dir> <subreads bam> <chrom> <start> <end> <fosmid name> <subreads_to_ref bam file>`

# ig - code for ig projects

[Fosmid related code](#fosmid-related-code)  

# Fosmid related code
### merge_seq_based_on_align_input.py 
```
python python/merge_seq_based_on_align_input.py <alignment_file> <fasta file with fosmid seq>
```
This is the code that was used to merge fosmid sequences and fosmids sequences with contigs from capture data. It takes as input the alignment file and a fasta file with sequences to merge. The alignmentfile is the file created in excel using the blast output. ## To do: Add example

# BAC related code
## assemble_barcoded_clones_without_backbone.sh
```
sh bash/assemble_barcoded_clones_without_backbone.sh \
  <outdir> \
  <input bam> \
  <genome size> \
  <backbone fasta>

Example pipeline:
bax1=XCWUL_20180912_RS42150_PL100105694-1_A01-1_lbc81*.1.bax.h5
bax2=XCWUL_20180912_RS42150_PL100105694-1_A01-1_lbc81*.2.bax.h5
bax3=XCWUL_20180912_RS42150_PL100105694-1_A01-1_lbc81*.3.bax.h5
bax2bam ${bax1} ${bax2} ${bax3} -o assemble/A1                                                                                                     
sh bash/assemble_barcoded_clones_without_backbone.sh \
  assemble \
  assemble/A1.subreads.bam \
  200000 \
  cloning_vector_pTARBAC2_1.fasta

```
This is the code that is used to assemble BACs.

