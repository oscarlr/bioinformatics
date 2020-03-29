#!/bin/env python
import sys
import pysam

tiling_path = sys.argv[1]
alignment = sys.argv[2]

def read_path(path):
    fosmids = {}
    with open(path, 'r') as pathfh:
        for line in pathfh:
            line = line.rstrip().split('\t')
            chrom = line[0]
            start = int(line[1])
            end = int(line[2])
            fosmid_name = line[3]
            fosmids[fosmid_name] = [chrom,start,end]
    return fosmids

def is_overlap(a, b):
    overlap =  max(0, min(a[1], b[1]) - max(a[0], b[0]))
    if overlap == 0:
        return False
    return True

# def overlapping_fosmids(path):
#     fosmids_that_overlap = []
#     fosmids = read_path(path)
#     for fosmid_1 in fosmids:
#         for fosmid_2 in fosmids:
#             if fosmid_1 == fosmid_2:
#                 continue
#             if is_overlap(fosmids[fosmid_1],fosmids[fosmid_2]):
#                 if sorted([fosmid_1,fosmid_2]) in fosmids_that_overlap:
#                     continue
#                 fosmids_that_overlap.append(sorted([fosmid_1,fosmid_2]))                                            
#     return fosmids_that_overlap

def read_completely_span(read,start,end):
    flank = 500
    if read.reference_start >= (start - flank):
        if read.reference_end <= (end + flank):
            return True
    return False

def get_number_of_reads(alignment,coord):
    samfile = pysam.AlignmentFile(alignment)
    number_of_reads = 0
    number_overlapping = 0
    for read in samfile.fetch(coord[0],coord[1],coord[2]):
        if not use_read(read):
            continue
        if read_completely_span(read,coord[1],coord[2]):
            number_of_reads += 1
        else:
            if is_overlap([read.reference_start,read.reference_end],[coord[1],coord[2]]):
                number_overlapping += 1
    return number_of_reads, number_overlapping

def get_reads_per_fosmid(alignment,path):
    fosmids = read_path(path)
    reads_per_fosmid = {}
    for fosmid in fosmids:
        number_of_reads,number_overlapping = get_number_of_reads(alignment,fosmids[fosmid])
        reads_per_fosmid[fosmid] = (number_of_reads,number_overlapping)
    return reads_per_fosmid

def use_read(read):
    if read.is_secondary:
        return False
    if read.is_supplementary:
        return False
    if read.is_duplicate:
        return False
    return True

def get_reads_outside_coords(alignment,path):
    fosmids = read_path(path)
    samfile = pysam.AlignmentFile(alignment)
    read_outside_coords = 0
    total_reads = 0
    for read in samfile:
        if not use_read(read):
            continue
        total_reads += 1
        read_span = False
        for fosmid in fosmids:
            if read_completely_span(read,fosmids[fosmid][1],fosmids[fosmid][2]):                
                read_span = True
                break
        if read_span:
            continue
        read_outside_coords += 1
    return (read_outside_coords,total_reads)

#fosmids = read_path(path)

#fosmids_that_overlap = overlapping_fosmids(tiling_path)
reads_per_fosmid = get_reads_per_fosmid(alignment,tiling_path)
read_outside_coords,total_reads = get_reads_outside_coords(alignment,tiling_path)

out_fosmids = []
for fosmids in fosmids_that_overlap:
    for fosmid in fosmids:
        if fosmid not in out_fosmids:            
            print "overlap\t%s" % fosmid
            out_fosmids.append(fosmid)

for fosmid in reads_per_fosmid:
    print "Number_of_reads\t%s\t%s\t%s" % (fosmid,reads_per_fosmid[fosmid][0],reads_per_fosmid[fosmid][1])

print "Reads_outside_coords\t%s\t%s\t%s" % (read_outside_coords,total_reads,float(read_outside_coords)/total_reads)
        
