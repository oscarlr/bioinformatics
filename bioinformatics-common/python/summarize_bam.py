#!/bin/env python
import sys
from scipy import stats
from scipy import ndimage
import numpy
import pysam

numbers = []
samfile = pysam.AlignmentFile(sys.argv[1],'rb',check_sq=False)
for read in samfile:
    numbers.append(int(read.query_length))
        
describe = stats.describe(numbers)
print "\t".join(["nobs","min","max","mean","std","median","coverage"])

out = [describe.nobs,
       describe.minmax[0],
       describe.minmax[1],
       describe.mean,
       ndimage.standard_deviation(numpy.array(numbers)),
       ndimage.median(numbers),
       ndimage.sum(numpy.array(numbers))/float(3100000000)]

print "\t".join(map(str,out))
