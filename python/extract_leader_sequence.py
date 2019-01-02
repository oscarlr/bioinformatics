#!/bin/env python
import sys
import pysam
from Bio.Seq import Seq
import re

fastafile = sys.argv[1]
sequencelistfn = sys.argv[2]

fasta = pysam.FastaFile(fastafile)
upstream_amount = 200

## Annotations can be in:
## ATG --> [GTX ---> XAG] --> [TAA,TGA,TAG]
## GT/AG in intron
## 

## Algorithm
## 1. Search for ATG upstream start of V gene segment
## 2. For each ATG search for all GT
## 3. Select searches that are in frame, ie count by 3s until you get to GT,
##    if no GT is found, discard


def read_in_genes(genesfn):
    genes = {}
    with open(genesfn,'r') as genesfh:
        for line in genesfh:
            line = line.rstrip().split('\t')
            gene_name = line[0]
            contig_name = line[1]
            start = int(line[2])
            end = int(line[3])
            genes[gene_name] = {}
            genes[gene_name]["location"] = [contig_name,start,end]
            genes[gene_name]["leader_sequences"] = None
    return genes

def extract_sequence(fasta,ref,start,end):
    seq = fasta.fetch(ref,start,end)
    return seq

def get_leader_seqs(seq):
    leader_seqs = []
    atg_start = [atg.start() for atg in re.finditer("ATG",seq)]
    for atg in atg_start:
        for i in range(atg,len(seq),3):
            frame = seq[i:i+3]
            if frame[0:2] == "GT":
                leader_seqs.append(seq[atg:i])
    return leader_seqs

genes = read_in_genes(sequencelistfn)
for gene in genes:
    ref = genes[gene]["location"][0]
    start = genes[gene]["location"][1] - upstream_amount
    end = genes[gene]["location"][2]
    gene_seq = extract_sequence(fasta,ref,start,end)
    genes[gene_name]["leader_sequences"] = get_leader_seqs(gene_seq)

                                
