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
