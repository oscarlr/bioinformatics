#!/bin/env python
import sys
import pysam

coordsfn = sys.argv[1]
bamfile = sys.argv[2]
coords = []

def get_read_coords(read,start,end):
    read_start = 0    
    read_end = None
    for query_pos, ref_pos in read.get_aligned_pairs():
        if query_pos == None:
            continue
        if ref_pos == None:
            continue
        if ref_pos < start:
            read_start = query_pos
        if ref_pos <= end:
            read_end = query_pos
        else:
            break
    assert read_end != None
    return (read_start,read_end)

with open(coordsfn,'r') as fh:
    for line in fh:
        line = line.rstrip().split('\t')
        chrom = line[0]
        start = int(line[1])
        end = int(line[2])
        #name = line[3]
        coords.append([chrom,start,end])#,name])

samfile = pysam.AlignmentFile(bamfile)
for chrom, start, end in coords:
    for read in samfile.fetch(chrom,start,end):
        if read.is_secondary:
            continue
        if read.is_supplementary:
            continue
        if read.is_unmapped:
            continue
        if read.reference_end < start:
            continue
        if read.reference_start > end:
            continue
        read_start, read_end = get_read_coords(read,start,end)
        read_sequence = read.query_sequence[read_start:read_end]
        print ">%s" % read.query_name
        print read_sequence

