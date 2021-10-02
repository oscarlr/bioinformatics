#!/bin/env python
import sys
import pysam

inbamfile = sys.argv[1]
outbamfile = sys.argv[2]
read_group = sys.argv[3]

insamfile = pysam.AlignmentFile(inbamfile,'rb')

header = dict(insamfile.header)
header['RG'] = [{'ID': read_group, 'SM': read_group}]

outsamfile = pysam.AlignmentFile(outbamfile,'wb',header=header)

def add_read_group(read,read_group):
    read_tags = read.get_tags()
    tags_to_add = []
    for tag in read_tags:
        if tag[0] != "RG":
            tags_to_add.append(tag)
    haptag = ("RG", read_group, "Z")
    tags_to_add.append(haptag)
    read.set_tags(tags_to_add)
    return read
    
for read in insamfile.fetch():
    if read.is_secondary:
        continue
    if read.is_unmapped:
        continue
    if read.is_supplementary:
        continue
    out_read = add_read_group(read,read_group)
    outsamfile.write(out_read)

insamfile.close()
outsamfile.close()
