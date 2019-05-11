#!/bin/python
import pysam
import sys

fasta = pysam.FastaFile(sys.argv[1])
entry_names = fasta.references
for entry_name in entry_names:
    sequence = fasta.fetch(reference=entry_name)
    sequence = set(list(sequence))
    non_fasta_letters = sequence.difference(set(["A","T","C","G","N","a","t","c","g","n"]))

print non_fasta_letters
    
