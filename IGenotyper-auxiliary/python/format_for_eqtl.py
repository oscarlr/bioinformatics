#!/bin/env python
import sys
import pysam

def read_samples(fofn):
    samples = {}
    with open(fofn,'r') as fh:
        for line in fh:
            line = line.rstrip().split('\t')
            sample_name = line[0]
            bam_file = line[1]
            vcf_file = line[2]
            samples[sample_name] = [bam_file,vcf_file]
    return samples

def get_non_deletion_pos(bam_file):
    """Returns genomic positions with aligned bases, eg no deletion bases
    Args:
        bam_file (str): BAM file path
    
    Returns:
        set: set of ints
    """
    pos = set()
    samfile = pysam.AlignmentFile(bam_file)
    for pileupcolumn in samfile.pileup("igh"):
        for pileupread in pileupcolumn.pileups:
            if pileupread.is_del:
                continue
            if pileupread.is_refskip:
                continue
            pos.add(pileupcolumn.pos)
    return pos

def load_snp_matrix(samples):
    snps = {}
    for sample in samples:
        bam_file, vcf_file = samples[sample]
        snps[sample] = {}
        non_deletion_pos = get_non_deletion_pos(bam_file)
        for pos in non_deletion_pos:
            snps[sample][pos] = "00"
    return snps

def get_all_positions(snp_matrix):
    all_positions = set()
    for sample in samples:
        all_positions.update(set(snp_matrix[sample].keys()))
    return all_positions

def load_nas(samples,snp_matrix):
    all_positions = get_all_positions(snp_matrix)
    for sample in samples:
        for pos in all_positions:
            if pos in sample:
                continue
            snp_matrix[sample][pos] = "NA"
    return snp_matrix
    
def get_genotype_pos(vcf_file):
    het_pos = set()
    hom_pos = set()
    with open(vcf_file,'r') as fh:
        for snpline in fh:
            if snpline[0] == "#":
                continue
            if "read_support=No" in snpline:
                continue
            if "sv_filter=Yes" in snpline:
                continue
            snpline = snpline.rstrip().split('\t')
            pos = int(snpline[1])
            genotype = snpline[9].split(":")[0]
            if genotype == "0/1":
                het_pos.add(pos)
            if genotype == "1/1":
                hom_pos.add(pos)
    return (het_pos,hom_pos)
                
def add_snps_to_matrix(samples,snp_matrix):
    for sample in samples:
        bam_file, vcf_file = samples[sample]
        het_pos, hom_pos = get_genotype_pos(vcf_file)
        for pos in het_pos:
            snp_matrix[sample][pos] = "10"
        for pos in hom_pos:
            snp_matrix[sample][pos] = "11"
    return snp_matrix

def write_matrix_to_file(snp_matrix,outfn):
    with open(outfn,'w') as fh:        
        out = ["snp_start"]
        for sample in snp_matrix:
            out.append(sample)
        fh.write("\t".join(out))
    all_positions = get_all_positions(snp_matrix)
    for pos in sorted(all_positions):
        out = [pos]
        for sample in snp_matrix:
            out.append(snp_matrix[sample][pos])
        fh.write("\t".join(map(str,out)))

def main():
    fofn = sys.argv[1]
    outfn  =sys.argv[2]
    samples = read_samples(fofn)
    snp_matrix = load_snp_matrix(samples)
    snp_matrix = load_nas(samples,snp_matrix)
    snp_matrix = add_snps_to_matrix(samples,snp_matrix)
    write_matrix_to_file(snp_matrix,outfn)
