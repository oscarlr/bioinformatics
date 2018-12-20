#!/bin/env python
import sys
import pysam

samfile = pysam.AlignmentFile(sys.argv[1])

for read in samfile:
    ref = samfile.get_reference_name(read.reference_id)
    start = read.reference_start
    end = read.reference_end
    out = [ref,start,end]
    print "%s" % "\t".join(map(str,out))
