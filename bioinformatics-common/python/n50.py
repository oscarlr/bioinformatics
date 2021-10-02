#!/bin/env python
import sys

faifn = sys.argv[1]
genome_size = float(sys.argv[2])

lengths = []
with open(faifn,'r') as fh:
    for line in fh:
        line = line.rstrip().split('\t')
        lengths.append(int(line[1]))

lengths = sorted(lengths, reverse = True)

cum_sum = 0
for l in lengths:
    cum_sum = cum_sum + l
    if cum_sum > (genome_size/2):
        assert l in lengths
        print(l)
        sys.exit(0)

print(None)
