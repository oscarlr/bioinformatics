#!/bin/env python
import sys
import vcf
import pysam

def read_phased_snps(vcffile,sample_name):
    phased_snps = {}
    vcf_reader = vcf.Reader(open(vcffile, 'r'))
    for record in vcf_reader:
        if record.var_subtype == "deletion":
            continue
        if record.var_subtype == "insertion":
            continue
        phased_snps.setdefault(record.CHROM, {})
        alleles  = record.genotype(sample_name)['GT'].split("|")
        if len(alleles) == 2 and alleles[0] != alleles[1]:
            alleles = map(int, alleles)
            if max(alleles) != 1:
                continue
            allele_bases = [record.REF,record.ALT[0]]
            sample_bases = map(lambda x: allele_bases[x], alleles)
            phased_snps[record.CHROM][record.POS - 1] = sample_bases
    return phased_snps

def calculate_no_prob(basetups,hap_dict):    
    b1_bases = 0.0
    b0_bases = 0.0
    call_per_pos = {}
    for basetup in basetups:
        rpos, b, qual = basetup
        b0, b1 = hap_dict[rpos]
        if b == b1 or b == b0:
            if b == b0:
                b0_bases += 1.0
                call_per_pos[rpos] = "b0"                
            else:
                b1_bases += 1.0                
                call_per_pos[rpos] = "b1"
    print "%s\t%s" % (b0_bases,b1_bases)

def phase_read(read,phased_snps,chrom):
    base_tuples = []
    if chrom in phased_snps:
        for read_pos, ref_pos in read.get_aligned_pairs():
            if ref_pos is None:
                continue
            if ref_pos in phased_snps[chrom]:
                if read_pos is not None:
                    qual = read.query_qualities[read_pos]
                    base = read.query_sequence[read_pos]
                    base_tuples.append((ref_pos,base,qual))
    if len(base_tuples) > 0:
        calculate_no_prob(base_tuples,phased_snps[chrom])

agreements_per_snp = {}
phased_snps = read_phased_snps(sys.argv[1],sys.argv[2])
unphased_bam = pysam.AlignmentFile(sys.argv[3],'rb')
for read in unphased_bam.fetch("igh",1,1179808):
    if read.is_secondary:
        continue
    if read.is_unmapped:
        continue
    if read.is_supplementary:
        continue
    if read.is_secondary:
        continue
    phase_read(read,phased_snps,unphased_bam.get_reference_name(read.reference_id))


