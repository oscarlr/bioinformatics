#!/bin/env python
import sys
import pysam
import itertools

bamfile = sys.argv[1]

def get_positions_with_bases(bamfile):
    positions = set()
    samfile = pysam.AlignmentFile(bamfile,"rb")
    for pileupcolumn in samfile.pileup("igh"):
        for pileupread in pileupcolumn.pileups:
            if not pileupread.is_del and not pileupread.is_refskip:
                positions.add(pileupcolumn.pos)
    return list(positions)

def turn_to_intervals(positions):
    for a, b in itertools.groupby(enumerate(positions), lambda pair: pair[1] - pair[0]):
        b = list(b)
        yield b[0][1], b[-1][1]
        

positions = get_positions_with_bases(bamfile)
for interval in  list(turn_to_intervals(positions)):
    out = ["igh",interval[0],interval[1]]
    print "\t".join(map(str,out))
