#!/bin/env python
import sys

vcffile = sys.argv[1]

def is_error(son,parent1,parent2):
    error = False
    if son == "0/1":
        if parent1 == "." and parent2 == ".":
            error = True
        if parent1 == "1/1" and parent2 == "1/1":
            error = True
    if son == "1/1":
        if parent1 == "." or parent2 == ".":
            error = True
    if son == ".":
        if parent1 == "1/1" or parent2 == "1/1":
            error = True
    return error

with open(vcffile,'r') as vcffh:
    for line in vcffh:
        if line.startswith("#"):
            continue
        line = line.rstrip().split('\t')
        son = line[9].split(":")[0]
        parent1 = line[10].split(":")[0]
        parent2 = line[11].split(":")[0]
        if is_error(son,parent1,parent2):
            print "\t".join(line)
