#!/bin/env python
from itertools import groupby, count
import pysam
import sys

PADDING = 1000
#start = int(sys.argv[2])
#end = int(sys.argv[3])

# def as_range(iterable): # not sure how to do this part elegantly
#     l = list(iterable)
#     if len(l) > 1:
#         return '{0}-{1}'.format(l[0], l[-1])
#     else:
#         return '{0}'.format(l[0])

regions = {}
with open(sys.argv[2]) as fh:
    for line in fh:
        line = line.rstrip().split('\t')
        ref = line[0]
        start = int(line[1])
        end = int(line[2])
        if ref not in regions:
            regions[ref] = []
        regions[ref].append((start,end))

numberlist = []
quals = []
with pysam.FastxFile(sys.argv[1]) as fh:
    for entry in fh:
        qual_array = entry.get_quality_array()
        for i,qual in enumerate(qual_array):
            add = False
            for start,end in regions[entry.name]:
                if i > start and i < end:
                    add = True
            if add:
                if qual < 40:
                    numberlist.append(i)
                    quals.append(qual)

#print quals
print len(numberlist)                
#print numberlist
#print ','.join(as_range(g) for _, g in groupby(numberlist, key=lambda n, c=count(): n-next(c)))
