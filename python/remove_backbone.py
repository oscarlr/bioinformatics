#!/bin/env python
import sys
import pysam

if __name__ == '__main__':
    samfile = pysam.AlignmentFile(sys.argv[1])
    for read in samfile:
        if read.is_unmapped:
            continue
        if read.is_secondary:
            continue
        if read.is_supplementary:
            continue

        left = read.query_sequence[:read.query_alignment_start]
        right = read.query_sequence[read.query_alignment_end:]
        outname = read.query_name
        if len(left) > 0:
            print ">%s_L\n%s" % (outname,left)
        if len(right) > 0:
            print ">%s_R\n%s" % (outname,right)

