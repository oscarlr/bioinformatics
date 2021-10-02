#!/bin/env python
import sys
from Bio import SeqIO

db = sys.argv[1]
species = sys.argv[2]
region = sys.argv[3]

with open(db,"r") as fh:
    for line in fh:
        line = line.rstrip()
        if ">" in line:            
            if species in line:                
                if region in line:                    
                    print line.replace(" ","")
                    print_seq = True
                else:
                    print_seq = False    
            else:
                print_seq = False
        else:
            if print_seq:
                print line
