#!/bin/env python
import sys

son_vcffile = sys.argv[1]
parent1_vcffile = sys.argv[2]
parent2_vcffile = sys.argv[3]

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

def read_genotype(vcffile):
    genotypes = {}
    with open(vcffile,'r') as vcffh:
        for line in vcffh:
            if line.startswith("#"):
                continue
            line = line.rstrip().split('\t')
            pos = line[1]
            genotype = line[9].split(":")[0]
            genotypes[pos] = genotype
    return genotypes

son_genotype = read_genotype(son_vcffile)
parent1_genotype = read_genotype(parent1_vcffile)
parent2_genotype = read_genotype(parent2_vcffile)

number_of_het = 0
number_of_het_errors = 0
number_of_hom = 0
number_of_hom_errors = 0
for pos in son_genotype:
    if pos not in parent1_genotype:
        continue
    if pos not in parent2_genotype:
        continue    
    if son_genotype[pos] == "0/1":
        number_of_het += 1.0        
    if son_genotype[pos] == "1/1":
        number_of_hom += 1.0      
    if is_error(son_genotype[pos],parent1_genotype[pos],parent2_genotype[pos]):
        print "error: %s" % pos
        if son_genotype[pos] == "0/1":
            number_of_het_errors += 1.0
        if son_genotype[pos] == "1/1":
            number_of_hom_errors += 1.0

print number_of_het
print number_of_hom
print number_of_het_errors/number_of_het
print number_of_hom_errors/number_of_hom

