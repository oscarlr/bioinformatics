#!/bin/env python
import sys

true_vcf = sys.argv[1]
test_vcf = sys.argv[2]
length = int(sys.argv[3])

def read_snp_pos(vcf):
    pos = set()
    with open(vcf,'r') as fh:
        for line in fh:
            if "#" in line:
                continue
            if "read_support=Yes" not in line:
                continue
            if "sv_filter=No" not in line:
                continue
            line = line.rstrip().split('\t')
            pos.add(line[1])
    return pos

true_snps = read_snp_pos(true_vcf)
test_snps = read_snp_pos(test_vcf)


true_positive = len(true_snps.intersection(test_snps))
false_negative = len(true_snps.difference(test_snps))
false_positive = len(test_snps.difference(true_snps))
true_negative = length - true_positive - false_positive

tpr = float(true_positive)/(true_positive + false_negative)
fpr = float(false_positive)/(false_positive + true_negative)

output = [tpr,fpr,true_positive,false_negative,false_positive,true_negative]
print "\t".join(map(str,output))
