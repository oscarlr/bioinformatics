#!/bin/env python
import sys
import pysam

samfile = pysam.AlignmentFile(sys.argv[1])
chrom = sys.argv[2]
start = int(sys.argv[3])
end = int(sys.argv[4])

for read in samfile.fetch(chrom,start,end):
    if read.is_unmapped:
        continue
    if read.is_secondary:
        continue
    if read.is_supplementary:
        continue
    if read.reference_start > end:
        continue
    if read.reference_end < start:
        continue
    print ">%s\n%s" % (read.query_name,read.query_sequence)
    
