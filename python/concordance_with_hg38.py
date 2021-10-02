#!/bin/env python
import sys
import pysam

# hg38
chrom = "chr14"
start = 105868036
end = 106871962

ref = pysam.FastaFile("/sc/orga/projects/bashia02c/projects/igh/refs/GR38/chr14_gr38.fasta")

seq = 0
indels = 0
mismatches = 0
alignment_length = 0

samfile = pysam.AlignmentFile(sys.argv[1], "rb" )
positions = set()
for pileupcolumn in samfile.pileup(chrom, start,end):
    if pileupcolumn.pos not in positions:
        seq += 1
        positions.add(pileupcolumn.pos)
    ref_base = ref.fetch(chrom,pileupcolumn.pos,pileupcolumn.pos + 1).upper()
    for pileupread in pileupcolumn.pileups:
        alignment_length += 1
        if pileupread.is_del or pileupread.is_refskip:
            indels += 1
            continue
        read_base = pileupread.alignment.query_sequence[pileupread.query_position].upper()
        if read_base != ref_base:
            mismatches += 1

no_seq = (end-start) - seq
perc_seq = float(seq)/(end-start)
accuracy = (alignment_length - (mismatches + indels))/float(alignment_length)

print "seq\t%s" % seq
print "indels\t%s" % indels
print "mismatches\t%s" % mismatches 
print "alignment_length\t%s" % alignment_length 
print "no_seq\t%s" % no_seq
print "perc_seq\t%s" % perc_seq
print "accuracy\t%s" % accuracy
