#!/bin/env python
import sys
import pysam
from collections import namedtuple

Alignment = namedtuple(
    'seq1',
    'seq2',
    'alignment_length',
    'perc_match',
    'mismatches',
    'num_gap_openings',
    'gap_length',
    'seq1_start',
    'seq1_end',
    'seq1_length',
    'seq2_start',
    'seq2_end',
    'seq2_length',
    'skip',
    'use',
    'new_name')

def get_alignments(input_file):
    alignments = []
    with open(input_file,'r') as input_fh:
        header = input_fh.readline()
        for line in input_fh:
            line = line.rstrip().split('\t')
            alignment = Alignment._make(line)
            alignments.append(alignment)
    return alignments

def read_sequences(input_fasta):
    sequences = {}
    fastafile = pysam.FastaFile(input_fasta)
    references = fastafile.references
    for reference in references:
        sequence = fastafile.fetch(reference,0,fastafile.get_reference_length(reference))
        sequences[reference] = sequence
    return sequences
    

def merge_sequences(alignments,sequences):
    for alignment in alignments:
        if alignment.use == "no":
            continue
        print ">%s" % alignment.new_name
        if alignment.seq2 == ".":
            new_seq = sequences[alignment.seq1]
        else:
            if alignment.seq1_start > alignment.seq2_start:
                new_seq = sequences[alignment.seq1] + sequences[alignment.seq2][alignment.seq2_end:]
            else:
                new_seq = sequences[alignment.seq2][:alignment.seq2_start] + sequences[alignment.seq1]
        print "%s" % new_seq
                
alignments = get_alignments(sys.argv[1])
sequences = read_sequences(sys.argv[2])
merge_sequences(alignments,sequences)
