#!/bin/env python
import sys
import pysam
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO
import datetime

def is_overlapping(a, b):
    if a[0] != b[0]: # chrom not matching
        return False
    overlapping = False
    num_overlapping = max(0, min(a[2], b[2]) - max(a[1], b[1]))
    if num_overlapping > 0:
        overlapping = True
    return overlapping

def get_haplotype(read_name):
    return read_name.split("_")[1].split('=')[1]

def assembly_location(read_name):
    '''
    Return a list of the chrom, start and end 
    of regions that was assembled.
    '''
    read_origin = read_name.split("_")[0].split('=')[1]
    chrom = read_origin.split(":")[0]
    start = int(read_origin.split(":")[1].split("-")[0])
    end = int(read_origin.split(":")[1].split("-")[1])
    return [chrom,start,end]

def snps_per_hap(bamfile,reffn,filter_on_region=True):
    samfile = pysam.AlignmentFile(bamfile)
    ref = pysam.FastaFile(reffn)
    regions = {}
    for read in samfile:
        if read.is_unmapped:
            continue
        if read.is_supplementary:
            continue
        if read.is_secondary:
            continue        
        assembled_region = assembly_location(read.query_name)
        mapped_chrom = samfile.get_reference_name(read.reference_id)
        mapped_start = int(read.reference_start)
        mapped_end = int(read.reference_end)
        mapped_region = [mapped_chrom,mapped_start,mapped_end]
        if filter_on_region:
            if not is_overlapping(assembled_region,mapped_region):
                continue
        haplotype = get_haplotype(read.query_name)
        for read_pos, ref_pos in read.get_aligned_pairs():
            if read_pos == None:
                continue
            if ref_pos == None:
                continue
            if filter_on_region:
                if ref_pos < int(assembled_region[1]):
                    continue
                if ref_pos > int(assembled_region[2]):
                    continue
            ref_base = ref.fetch(mapped_chrom,ref_pos,ref_pos + 1).upper()
            read_base = read.query_sequence[read_pos].upper()
            if ref_base != read_base:
                read_qual = 8 #read.query_qualities[read_pos]
                if mapped_chrom not in regions:
                    regions[mapped_chrom] = {}
                if ref_pos not in regions[mapped_chrom]:
                    regions[mapped_chrom][ref_pos] = {}                
                regions[mapped_chrom][ref_pos][haplotype] = (ref_base,read_base,read_qual,read.query_name)
    return regions


def add_genotype_to_snps(haplotype_snps):
    output = []
    for chrom in haplotype_snps:
        for pos in haplotype_snps[chrom]:
            if "haploid" in haplotype_snps[chrom][pos]:
                ref,read_base,read_qual,read_name = haplotype_snps[chrom][pos]["haploid"]
                line = [chrom,int(pos),".",ref,read_base,"./.",read_qual,read_name]
                output.append(line)
            else:
                if "0" in haplotype_snps[chrom][pos]:
                    ref,read_base,read_qual,read_name = haplotype_snps[chrom][pos]["0"]
                    line = [chrom,int(pos),".",ref,read_base,"1/1",read_qual,read_name]
                    output.append(line)
                else:
                    if len(haplotype_snps[chrom][pos]) == 1:
                        hap = None
                        if "1" in haplotype_snps[chrom][pos]:
                            hap = "1"
                        if "2" in haplotype_snps[chrom][pos]:
                            hap = "2"
                        ref,read_base,read_qual,read_name = haplotype_snps[chrom][pos][hap]
                        line = [chrom,int(pos),".",ref,read_base,"0/1",read_qual,read_name]
                        output.append(line)
                    else:
                        hap_1_ref,hap_1_read_base,hap_1_read_qual,hap_1_read_name = haplotype_snps[chrom][pos]["1"]
                        hap_2_ref,hap_2_read_base,hap_2_read_qual,hap_2_read_name = haplotype_snps[chrom][pos]["2"]
                        if hap_1_read_base == hap_2_read_base:
                            genotype = "1/1"
                            alt_base = hap_1_read_base
                        else:
                            genotype = "1/2"
                            alt_base = "%s,%s" % (hap_1_read_base,hap_2_read_base)
                        read_names = "%s,%s" % (hap_1_read_name,hap_2_read_name)
                        alt_qual = "%s,%s" % (hap_1_read_qual,hap_2_read_qual)
                        line = [chrom,int(pos),".",hap_1_ref,alt_base,genotype,alt_qual,read_names]
                        output.append(line)
    return output

def vcf_header(sample_name):
    i = datetime.datetime.now()
    line = [ "##fileformat=VCFv4.2",
             "##fileDate=%s%s%s" % (i.year,i.month,i.day),
             "##source=IGenotyper",
             "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">",
             "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
             "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">",
             "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s" % sample_name]
    return line

def write_snps_output_to_vcf(genotyped_snps,sample_name,vcffile):
    output_lines = vcf_header(sample_name)
    genotyped_snps.sort(key=lambda x: x[1])
    with open(vcffile,"w") as vcf_fh:
        vcf_fh.write("%s\n" % "\n".join(output_lines))
        for genotyped_snp in genotyped_snps:
            genotype_with_qual = [":".join(map(str,genotyped_snp[5:-1]))]            
            if type(genotyped_snp[-2]) != int:
                average_quality = sum(map(int,genotyped_snp[-2].split(",")))/2
            else:
                average_quality = genotyped_snp[-2]
            append_to_line = [average_quality,"PASS"]
            outline = genotyped_snp[:5] + append_to_line + genotype_with_qual
            vcf_fh.write("%s\n" % "\t".join(map(str,outline)))

mapped_locus = sys.argv[1]
ref = sys.argv[2]
sample_name = sys.argv[3]
vcffile = sys.argv[4]

haplotype_snps = snps_per_hap(mapped_locus,ref)
genotyped_snps = add_genotype_to_snps(haplotype_snps)
write_snps_output_to_vcf(genotyped_snps,sample_name,vcffile)
