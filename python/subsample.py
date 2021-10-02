#!/bin/env python
import sys
import pysam
import random

bamfile = sys.argv[1]
outdir = sys.argv[2]

subsample_values = [1,2,3,4,5,6,7,8,9]

subsample_bams = {}

samfile = pysam.AlignmentFile(bamfile, 'rb',check_sq=False)
header = samfile.header.copy()

for value in subsample_values:
    bamfile = "%s/bam_downsample_%s.bam" % (outdir,str(value))
    subsample_bams[value] = pysam.Samfile(bamfile, 'wb', header=header)
    
for read in samfile:
    random_value = random.randint(0,101)
    for value in subsample_values:
        if random_value < value:
            subsample_bams[value].write(read)

for value in subsample_bams:
    subsample_bams[value].close()
    
    
