#!/bin/env python
import sys
import pysam

bamfn = sys.argv[1]
outbamfn = sys.argv[2]
groupstxt = sys.argv[3]

def read_in_groups(groupstxt):
    groups = {}
    with open(groupstxt,'r') as fh:
        for line in fh:
            line = line.rstrip().split('\t')
            group = line[0]
            fosmid = line[1]
            if group not in groups:
                groups[group] = []
            groups[group].append(fosmid)
    return groups

def create_header(unphased_bam,groups):
    read_groups = []
    for group in groups:
        read_groups.append({"ID": group})
    if "RG" in unphased_bam.header:
        for group in unphased_bam.header["RG"]:
            for item in group:
                if item != "ID":
                    for ingroup in read_groups:
                        ingroup[item] = group[item]
    phased_bam_header = unphased_bam.header.to_dict()
    phased_bam_header["RG"] = read_groups
    return phased_bam_header

def write_bam(bamfn,outbamfn,groups):
    unphased_bam = pysam.AlignmentFile(bamfn,'rb')  
    header = create_header(unphased_bam,groups)
    phased_bam = pysam.AlignmentFile(outbamfn,'wb',header=header)
    for read in unphased_bam.fetch(): 
        if read.is_secondary:
            continue
        if read.is_unmapped:
            continue
        if read.is_supplementary:
            continue
        if read.is_secondary:
            continue
        read_tags = read.get_tags()
        tags_to_add = []
        for tag in read_tags:
            if tag[0] != "RG":
                tags_to_add.append(tag)        
        for group in groups:
            if read.query_name in groups[group]:                
                grouptag = ("RG", str(group), "Z")
                tags_to_add.append(grouptag)
        read.set_tags(tags_to_add)
        phased_bam.write(read)
    phased_bam.close()
    unphased_bam.close()

groups = read_in_groups(groupstxt)
write_bam(bamfn,outbamfn,groups)
