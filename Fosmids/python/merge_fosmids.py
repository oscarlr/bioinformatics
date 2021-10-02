#!/bin/env python
import sys
import networkx 
from collections import namedtuple
from networkx.algorithms.components.connected import connected_components

# from ..command_line import *
# from ..common import *
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# def run_blast_on_contigs(locus_fasta,locus_to_locus_blast):    
#     run_blast(locus_fasta,locus_fasta,blast_output)

def not_overlapping(alignment):
    flank = 500
    if not (int(alignment.qstart) < flank and (int(alignment.slen) - int(alignment.send)) < flank):
        if not (int(alignment.sstart) < flank and (int(alignment.qlen) - int(alignment.qend)) < flank):
            return True
    return False
   
def completely_overlapping(alignment):
    q_start = int(alignment.qstart) == 1
    q_end = (int(alignment.qlen) - int(alignment.qend)) == 0
    s_start = int(alignment.sstart) == 1
    s_end = (int(alignment.slen) - int(alignment.send)) == 0
    if (q_start and q_end) or (s_start and s_end):
        return True
    return False

def filter_alignments(locus_to_locus_blast,length,errors):
    pairs = []
    alignments = []
    pident_min = 95.0            
    mismatch_max = errors
    length_min = length # 10,000
    columns = ["length","pident","nident","mismatch","gapopen","gaps","qseqid",
               "qstart","qend","qlen","sseqid","sstart","send","slen","sstrand"]
    fosmid_names = set()
    Alignment = namedtuple('Alignment',columns)
    with open(locus_to_locus_blast,'r') as fh:
        for line in fh:
            line = line.rstrip().split('\t')
            alignment = Alignment._make(line)            
            fosmid_names.add(alignment.qseqid)
            fosmid_names.add(alignment.sseqid)
            if int(alignment.length) < length_min:
                continue
            if alignment.qseqid == alignment.sseqid:
                continue
            if float(alignment.pident) < pident_min:
                continue
            if float(alignment.mismatch) > mismatch_max:
                continue
            if alignment.sstrand != "plus":
                continue            
            if not_overlapping(alignment):
                continue
            if int(alignment.sstart) > int(alignment.qstart):
                continue
            if sorted((alignment.qseqid,alignment.sseqid)) in pairs:
                continue
            # if completely_overlapping(alignment):
            #     continue
            pairs.append(sorted((alignment.qseqid,alignment.sseqid)))
            alignments.append(alignment)
    return alignments,fosmid_names

def to_graph(l):
    G = networkx.Graph()
    for part in l:
        # each sublist is a bunch of nodes
        G.add_nodes_from(part)
        # it also imlies a number of edges:
        G.add_edges_from(to_edges(part))
    return G

def to_edges(l):
    """ 
        treat `l` as a Graph and returns it's edges 
        to_edges(['a','b','c','d']) -> [(a,b), (b,c),(c,d)]
    """
    it = iter(l)
    last = next(it)

    for current in it:
        yield last, current
        last = current    

def group_alignments(alignments,fosmids_to_ignore,fosmid_names):
    fosmid_groupings = []
    for alignment in alignments:
        if alignment.qseqid in fosmids_to_ignore:
            continue
        if alignment.sseqid in fosmids_to_ignore:
            continue
        if alignment.qseqid in fosmid_names:
            fosmid_names.remove(alignment.qseqid)
        if alignment.sseqid in fosmid_names:
            fosmid_names.remove(alignment.sseqid)            
        fosmid_groupings.append([alignment.qseqid,alignment.sseqid])
    G = to_graph(fosmid_groupings)
    groupings = {}
    for i,group in enumerate(connected_components(G)):
        groupings[i] = group
    if len(groupings.keys()) != 0:
        max_group = max(groupings.keys())
    else:
        max_group = 0
    for i,fosmid in enumerate(fosmid_names):
        groupings[i + max_group + 1] = [fosmid]
    return groupings

def create_directed_graph(alignments,fosmids):
    G = networkx.DiGraph()
    for alignment in alignments:
        if alignment.qseqid not in fosmids:
            continue
        if alignment.sseqid not in fosmids:
            continue
        weight = int(alignment.length) + (int(alignment.slen) - int(alignment.length)) + (int(alignment.qlen) - int(alignment.length))
        G.add_edge(alignment.qseqid, alignment.sseqid, weight=weight)
    return G

def get_fosmid_coords(alignments,max_path):
    if len(max_path) == 1:
        #return [[max_path[0],1,-1]]
        return [[max_path[0],1,-1]]
    coords = []
    i = 0
    for start_fosmid, end_fosmid in zip(max_path,max_path[1:]):
        for alignment in alignments:
            if alignment.qseqid != start_fosmid:
                continue
            if alignment.sseqid != end_fosmid:
                continue
            merge_coords = [alignment.qseqid,1,int(alignment.qstart) - 1]
            #merge_coords = [alignment.qseqid,1,int(alignment.qlen)]
            coords.append(merge_coords)
            i += 1
            if i == len(max_path) - 1:                
                merge_coords = [alignment.sseqid,1,alignment.slen]
                #merge_coords = [alignment.sseqid,alignment.send,alignment.slen]
                coords.append(merge_coords)
    return coords
    
def group_size_2_and_encaps(group,alignments):
    if len(group) == 2:
        for alignment in alignments:
            if (alignment.qseqid == group[0] and alignment.sseqid == group[1]) or \
               (alignment.qseqid == group[1] and alignment.sseqid == group[0]):
                if completely_overlapping(alignment):
                    return True
    return False

def get_larger_contig(group,alignments):
    larger_contig = None
    for alignment in alignments:
        if (alignment.qseqid == group[0] and alignment.sseqid == group[1]) or \
           (alignment.qseqid == group[1] and alignment.sseqid == group[0]):
            if int(alignment.qlen) > int(alignment.slen):
                larger_contig = alignment.qseqid
            else:
                larger_contig = alignment.sseqid
            break
    assert larger_contig != None
    return larger_contig

def group_merging_fosmids(alignments,fosmids_to_ignore,fosmid_names):
    groupings = {}
    groups = group_alignments(alignments,fosmids_to_ignore,fosmid_names)
    for group in groups:        
        G = create_directed_graph(alignments,groups[group])        
        all_paths_weighted = dict(networkx.all_pairs_dijkstra_path_length(G))
        all_paths =  networkx.shortest_path(G)
        max_path = None
        num_fosmids_in_max_path = 0
        for start in all_paths_weighted:
            for end in all_paths_weighted[start]:
                length = all_paths_weighted[start][end]
                if length > num_fosmids_in_max_path:
                    num_fosmids_in_max_path = length
                    max_path = all_paths[start][end]
        if group_size_2_and_encaps(list(groups[group]),alignments):
            max_path = None
            groups[group] = [get_larger_contig(list(groups[group]),alignments)]
        if max_path == None:
            max_path = groups[group]
        fosmids_coordinates = get_fosmid_coords(alignments,max_path)
        groupings[group] = fosmids_coordinates
    return groupings

# def group_merging_fosmids(alignments,fosmids_to_ignore,fosmid_names):
#     groupings = {}
#     groups = group_alignments(alignments,fosmids_to_ignore,fosmid_names)
#     for group in groups:        
#         G = create_directed_graph(alignments,groups[group])        
#         all_paths =  networkx.shortest_path(G)
#         max_path = None
#         num_fosmids_in_max_path = 0
#         for start in all_paths:
#             for end in all_paths[start]:
#                 if len(all_paths[start][end]) > num_fosmids_in_max_path:
#                     num_fosmids_in_max_path = len(all_paths[start][end])
#                     max_path = all_paths[start][end]
#         if max_path == None:
#             max_path = groups[group]
#         fosmids_coordinates = get_fosmid_coords(alignments,max_path)
#         groupings[group] = fosmids_coordinates
#     return groupings

def read_merge_instructions(merge_alignments_instructions):
    columns = ["contig_group","igh_1","qseqid_start_1","qseqid_end_1","qseqid_hap_1",
               "igh_2","sseqid_start_1","sseqid_end_1","sseqid_hap_1",
               "length","pident","nident","mismatch","gapopen","gaps",
               "qseqid","qstart","qend","qlen",
               "sseqid","sstart","send","slen",
               "sstrand"]
    Alignment = namedtuple('Alignment',columns)
    groupings = {}
    with open(merge_alignments_instructions,'r') as fh:
        for line in fh:
            line = line.rstrip().split('\t')
            alignment = Alignment._make(line)            
            if alignment.contig_group not in groupings:
                groupings[alignment.contig_group] = []
                groupings[alignment.contig_group].append([alignment.qseqid,1,alignment.qlen,alignment.qseqid_hap_1])
                groupings[alignment.contig_group].append([alignment.sseqid,alignment.send,alignment.slen,alignment.sseqid_hap_1])
                continue
            assert alignment.qseqid == groupings[alignment.contig_group][-1][0]
            groupings[alignment.contig_group].append([alignment.sseqid,alignment.send,alignment.slen,alignment.sseqid_hap_1])
    return groupings

def merge_contigs(groupings,locus_fasta,merged_contigs,single_contigs_to_add):
    contigs = SeqIO.to_dict(SeqIO.parse(locus_fasta,"fasta"))
    contigs_to_output = []
    contigs_to_be_merged = [contigs[i][0] for i in contigs]    
    for group in groupings:
        hap=None
        starts=[]
        ends=[]
        sequences = []
        for contig_name,contig_start,contig_end,contig_hap in groupings[group]:            
            if hap == None or hap == "0":
                hap = contig_hap            
            locus_start,locus_end,c_hap = get_entry_location(contig_name)
            starts.append(int(locus_start))
            ends.append(int(locus_end))
            sequence = str(contigs[contig_name].seq[int(contig_start):int(contig_end)])
            sequences.append(sequence)
        sequence = "".join(sequences)
        sequence_name = "coord=igh:%s-%s_hap=%s_index=%s_total=merged_/0/0_0" % (min(starts),max(ends),hap,group)
        record = SeqRecord(Seq(sequence,"fasta"),id=sequence_name,name="",description="")
        contigs_to_output.append(record)
    with open(single_contigs_to_add,'r') as fh:
        for line in fh:
            contig = line.rstrip()
            contigs_to_output.append(contigs[contig])
    SeqIO.write(contigs_to_output,merged_contigs,"fasta")

def read_fosmids_to_ignore(fosmids_to_ignorefn):
    fosmids = []
    with open(fosmids_to_ignorefn,'r') as fh:
        for line in fh:
            line = line.rstrip().split('\t')
            fosmids.append(line[0])
    return fosmids

def save_alignments(alignments,alignments_file):
    with open(alignments_file,'w') as fh:
        for alignment in alignments:
            fh.write("%s\n" % "\t".join(map(str,alignment)))

def write_sequences(infosmids_fasta,outfosmids_fasta,fosmids_coords_to_mergefn):
    infasta = SeqIO.to_dict(SeqIO.parse(infosmids_fasta,"fasta"))
    outsequences = {}
    with open(fosmids_coords_to_mergefn,'r') as fh:
        for line in fh:
            line = line.rstrip().split('\t')
            group = line[0]
            fosmid = line[1]
            start = int(line[2])
            end = int(line[3])
            if group not in outsequences:
                outsequences[group] = []
            if start == 1 and end == -1:
                seq = str(infasta[fosmid].seq[0:])
            else:
                seq = str(infasta[fosmid].seq[start:end])            
            outsequences[group].append(seq)
    records = []
    for group in outsequences:
        sequence = "".join(outsequences[group])        
        record = SeqRecord(Seq(sequence,"fasta"),id="%s/0/0_8" % group,name="",description="")
        records.append(record)
    SeqIO.write(records,outfosmids_fasta,"fasta")

def get_possible_merges(fosmdis_to_fosmids_blast,blast_filtered,groupingsfn,fosmids_to_ignorefn,fosmids_coords_to_mergefn,infosmids_fasta,outfosmids_fasta,length,errors):
    alignments,fosmid_names = filter_alignments(fosmdis_to_fosmids_blast,length,errors)
    fosmids_to_ignore = read_fosmids_to_ignore(fosmids_to_ignorefn)
    groupings = group_alignments(alignments,fosmids_to_ignore,fosmid_names)
    groupings_with_fosmids_to_merge = group_merging_fosmids(alignments,fosmids_to_ignore,fosmid_names)
    save_alignments(alignments,blast_filtered)
    with open(groupingsfn,'w') as fh:
        for group in groupings:
            for fosmid_name in groupings[group]:
                fh.write("%s\t%s\n" % (group,fosmid_name))
    with open(fosmids_coords_to_mergefn,'w') as fh:
        for group in groupings_with_fosmids_to_merge:
            for fosmid_coord in groupings_with_fosmids_to_merge[group]: 
                output = [group] + fosmid_coord
                fh.write("%s\n" % "\t".join(map(str,output)))
    write_sequences(infosmids_fasta,outfosmids_fasta,fosmids_coords_to_mergefn)


fosmdis_to_fosmids_blast = sys.argv[1]
blast_filtered = sys.argv[2]
groupingsfn = sys.argv[3]
fosmids_to_ignorefn = sys.argv[4]
fosmids_coords_to_mergefn = sys.argv[5]
infosmids_fasta = sys.argv[6]
outfosmids_fasta = sys.argv[7]
length = int(sys.argv[8])
errors = int(sys.argv[9])
get_possible_merges(fosmdis_to_fosmids_blast,blast_filtered,groupingsfn,fosmids_to_ignorefn,fosmids_coords_to_mergefn,infosmids_fasta,outfosmids_fasta,length,errors)
