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
