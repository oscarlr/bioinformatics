#!/bin/env python
import sys
import pysam

samfile = pysam.AlignmentFile(sys.argv[1])

for read in samfile:
    ref = samfile.get_reference_name(read.reference_id)
    start = read.reference_start
    end = read.reference_end
    strand = "+"
    if read.is_reverse:
        strand = "-"
    out = [ref,start,end,strand]
    print "%s" % "\t".join(map(str,out))
