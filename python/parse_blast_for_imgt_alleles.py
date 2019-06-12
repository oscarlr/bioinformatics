#!/bin/env python
import sys
import pysam
from collections import namedtuple

columns =  ['seq1',
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
            'strand']

Alignment = namedtuple("Alignment",columns)

def get_alignments(input_file):
    alignments = []
    with open(input_file,'r') as input_fh:
        for line in input_fh:
            line = line.rstrip().split('\t')
            alignment = Alignment._make(line)
            alignments.append(alignment)   
    return alignments

def overlap(coordinate,list_of_coordinates):
    for coord_in_list in list_of_coordinates:
        if coordinate[0] != coord_in_list[0]:
            continue
        amount_overlap = max(0,min(int(coord_in_list[2]),int(coordinate[2])) - max(int(coordinate[1]),int(coord_in_list[1])))
        if amount_overlap > 0:
            return True
    return False

def seq_is_subset(seq,sequences):
    for seq_in_list in sequences:
        if seq == seq_in_list:
            continue
        if seq in seq_in_list:
            return True
    return False

def get_genes(fosmids):
    genes = set()
    for fosmid in fosmids:
        gene_hit = fosmids[fosmid][0].split("|")[1].split("*")[0]
    genes.add(gene_hit)
    return genes
        
alignments = get_alignments(sys.argv[1])
fasta = pysam.FastaFile(sys.argv[2])

matches = {}
perfect_coordinates = []
for alignment in alignments:
    key = (alignment.seq1,alignment.seq1_start,alignment.seq1_end)
    if float(alignment.perc_match) == 100.0:
        perfect_coordinates.append(key)
    if key not in matches:
        matches[key] = [alignment.seq2,float(alignment.perc_match),alignment.strand,alignment.alignment_length]
    else:
        if float(alignment.perc_match) > matches[key][1]:
            matches[key] = [alignment.seq2,float(alignment.perc_match),alignment.strand,alignment.alignment_length]

sequences = {}
for match in matches:
    if matches[match][1] == 100.0:
        continue
    if overlap(match,perfect_coordinates):
        continue
    sequence = fasta.fetch(list(match)[0],int(list(match)[1]),int(list(match)[2])).upper()
    assert (len(sequence) + 1) == int(matches[match][3])
    if sequence not in sequences:
        sequences[sequence] = {}
    gene_hit = matches[match][0]
    sequences[sequence][list(match)[0]] = [gene_hit,list(match)[1],list(match)[2]]

header = ["genes","seq","fosmid_to_output","fosmid_to_output_start","fosmid_to_output_end","number_of_fosmids","fosmids","fosmid_sequence"]
print "\t".join(map(str,header))
for seq in sequences:
    if seq_is_subset(seq,sequences.keys()):
        continue
    genes = ",".join(list(get_genes(sequences[seq])))
    fosmid_to_output = sequences[seq].keys()[0]
    fosmid_to_output_start = sequences[seq][fosmid_to_output][1]
    fosmid_to_output_end = sequences[seq][fosmid_to_output][2]
    fosmids = ",".join(sequences[seq].keys())
    number_of_fosmids = len(sequences[seq].keys())
    fosmid_sequence = fasta.fetch(fosmid_to_output)
    out = [genes,seq,fosmid_to_output,fosmid_to_output_start,fosmid_to_output_end,number_of_fosmids,fosmids,fosmid_sequence]
    print "\t".join(map(str,out))
