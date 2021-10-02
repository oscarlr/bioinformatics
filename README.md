# bioinformatics
## python scripts

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

