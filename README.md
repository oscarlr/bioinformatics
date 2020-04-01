# bioinformatics
## python scripts

### Get sequence (fasta file) from bam file for each bed entry
`python/extract_sequence_from_bam.py <bed file> <bam file> > <fasta file>`

### Get tab delimited file with number of reads per bed entry file
`python/check_mapping_in_coords.py <bed file> <bam file> > <tab-delim file>`

### Get regions (bed file) with coverage
`python/get_pos_with_bases.py <bam file> > <bed file>`
