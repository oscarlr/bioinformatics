#!/bin/env python
import os
import sys
import time
import requests
from Bio import SeqIO
from io import StringIO
from vquest import util,request
from collections import namedtuple
from zipfile import ZipFile
import argparse
#from vquest import * # https://github.com/ressy/vquest

URL = "https://www.imgt.org/IMGT_vquest/analysis"

CONFIG = {'outputType': 'html',
          'inputType': 'inline',
          'nbNtPerLine': 60,
          'resultType': 'excel',
          'dv_V_GENEalignment': True,
          'dv_D_GENEalignment': False,
          'dv_J_GENEalignment': True,
          'dv_IMGTjctaResults': True,
          'dv_eligibleD_GENE': False,
          'dv_JUNCTIONseq': True,
          'dv_V_REGIONalignment': True,
          'dv_V_REGIONtranlation': True,
          'dv_V_REGIONprotdisplay': True,
          'dv_V_REGIONmuttable': True,
          'dv_V_REGIONmutstats': True,
          'dv_V_REGIONhotspots': True,
          'dv_IMGTgappedVDJseq': True,
          'dv_IMGTAutomat': True,
          'dv_IMGTCollierdePerles': 0,
          'sv_V_GENEalignment': True,
          'sv_V_REGIONalignment': True,
          'sv_V_REGIONtranslation': True,
          'sv_V_REGIONprotdisplay': True,
          'sv_V_REGIONprotdisplay2': True,
          'sv_V_REGIONprotdisplay3': True,
          'sv_V_REGIONfrequentAA': True,
          'sv_IMGTjctaResults': True,
          'xv_outputtype': 3,
          'xv_summary': True,
          'xv_IMGTgappedNt': True,
          'xv_ntseq': True,
          'xv_IMGTgappedAA': True,
          'xv_AAseq': True,
          'xv_JUNCTION': True,
          'xv_V_REGIONmuttable': True,
          'xv_V_REGIONmutstatsNt': True,
          'xv_V_REGIONmutstatsAA': True,
          'xv_V_REGIONhotspots': True,
          'xv_parameters': True,
          'IMGTrefdirSet': 1,
          'IMGTrefdirAlleles': True,
          'V_REGIONsearchIndel': True,
          'nbD_GENE': -1,
          'nbVmut': -1,
          'nbDmut': -1,
          'nbJmut': -1,
          'nb5V_REGIONignoredNt': 0,
          'nb3V_REGIONaddedNt': 0,
          'scfv': False,
          'cllSubsetSearch': False,
          'species': None, # user-specified
          'receptorOrLocusType': None, # user-specified
          'fileSequences': None, # user-specified
          'sequences': None # user-specified
          }

CHUNK_SIZE = 50

DELAY = 1

def parse_gene_type(name):
    # feat=TRBV2_hap=2_pos=chr7:142301134-142301432_i=0
    gene_type = name.split('_')[0]
    gene_type = gene_type.split('=')[1]
    gene_type = gene_type[3]
    return gene_type

def parse_gene(name):
    # feat=TRBV2_hap=2_pos=chr7:142301134-142301432_i=0
    gene = name.split('_')[0]
    gene = gene.split('=')[1]
    return gene

def parse_hap(name):
    hap = name.split('_')[1]
    hap = hap.split('=')[1]
    return hap

def get_element(imgt,element,value):
    if element in imgt:
        value = imgt[element]
    return value

def get_gene_element(imgt,element,gene):
    value = None
    element = "%s_%s" % (gene,element)
    if element in imgt:
        value = imgt[element]
    return value
        
class ValidateAllele:
    def __init__(self,name,seq):
        self.sequence_id = name
        self.seq = seq
        self.gene = parse_gene(self.sequence_id)
        self.gene_type = parse_gene_type(self.sequence_id)
        self.hap = parse_hap(self.sequence_id) 
        
        # IMGT
        self.sequence = None
        self.productive = None
        self.call = None
        self.identity = None
        self.sequence_alignment = None
        self.germline_alignment = None
        self.full_length = None
        
        # Additional support
        self.ccs_support = 0
        self.in_other_samples = None

    def add_imgt_data(self,imgt):
        self.sequence = get_element(imgt,"sequence",self.sequence)
        self.productive = get_element(imgt,"productive",self.productive)
        self.call = get_gene_element(imgt,"call",self.gene_type.lower())                                     
        self.identity = get_gene_element(imgt,"identity",self.gene_type.lower())
        self.sequence_alignment = get_gene_element(imgt,"sequence_alignment",self.gene_type.lower())
        self.germline_alignment = get_gene_element(imgt,"germline_alignment",self.gene_type.lower())
        if self.sequence_alignment not in [None,""]:
            if "." == self.sequence_alignment[0]:
                self.full_length = False
            elif "." == self.sequence_alignment[-1]:
                self.full_length = False
            else:
                self.full_length = True
                
    def add_ccs_support(self,ccs_counts):
        if self.gene in ccs_counts:
            if self.hap in ccs_counts[self.gene]:
                if self.seq in ccs_counts[self.gene][self.hap]:
                    self.ccs_support = ccs_counts[self.gene][self.hap][self.seq]

    def found_in_other_samples(self,other_samples):
        if self.seq in other_samples:
            self.in_other_samples = other_samples[self.seq]
        else:
            self.in_other_samples = 0
            
    def __repr__(self):
        out = [self.sequence_id,
               self.gene_type,
               self.gene,
               self.hap,
               self.call,
               self.identity,
               self.productive,
               self.ccs_support,
               self.in_other_samples,
               self.full_length,
               self.sequence,
               self.sequence_alignment,
               self.germline_alignment,
               self.seq]
        out = "\t".join(map(str,out))
        return out

    def __str__(self):
        return self.__repr__()
        

def ccs_read_count(fn):
    ccs_alleles = {}
    for seq_record in SeqIO.parse(fn, "fasta"):
        gene = parse_gene(seq_record.id)
        hap = parse_hap(seq_record.id)
        seq = str(seq_record.seq)
        if gene not in ccs_alleles:
            ccs_alleles[gene] = {}
        if hap not in ccs_alleles[gene]:
            ccs_alleles[gene][hap] = {}
        if seq not in ccs_alleles[gene][hap]:
            ccs_alleles[gene][hap][seq] = 0
        ccs_alleles[gene][hap][seq] += 1
    return ccs_alleles

def read_other_samples(fofn):
    seqs = {}
    with open(fofn,'r') as fofh:
        for line in fofh:
            fn = line.rstrip()
            for seq_record in SeqIO.parse(fn, "fasta"):
                if str(seq_record.seq) not in seqs:
                    seqs[str(seq_record.seq)] = 0
                seqs[str(seq_record.seq)] += 1
    return seqs

def set_run_config(species,receptorOrLocusType,fileSequences):
    run_config = CONFIG
    run_config['species'] = species
    run_config['receptorOrLocusType'] = receptorOrLocusType
    run_config['fileSequences'] = fileSequences
    return run_config

def set_chunk_config(run_config,chunk):
    config_chunk = run_config.copy()
    out_handle = StringIO()
    SeqIO.write(chunk, out_handle, "fasta")
    config_chunk["sequences"] = out_handle.getvalue()    
    return config_chunk

def run_vquest(run_config,assemble_alleles,imgt_outputfn):
    added_header = False
    count = 0
    if not os.path.isfile(imgt_outputfn):
        outfn = open(imgt_outputfn,'w')
        for i,chunk in enumerate(util.chunker(assemble_alleles, CHUNK_SIZE)):
            chunk_config = set_chunk_config(run_config,chunk)    
            response = requests.post(URL, data = chunk_config)
            with open("_tmp.zip", 'wb') as f:
                f.write(response.content)
            with ZipFile("_tmp.zip", 'r') as zip:
                data = zip.read("vquest_airr.tsv")
            data = data.decode().split('\n')
            for j,line in enumerate(data):                
                if j == 0:
                    if added_header:
                        continue
                    else:
                        added_header = True
                outfn.write("%s\n" % line)
            time.sleep(DELAY)
        outfn.close()
        
def parse_vquest(imgt_outputfn):
    header = None
    vquest_output = {}
    with open(imgt_outputfn,'r') as fh:
        for i,line in enumerate(fh):
            line = line.rstrip().split('\t')
            if len(line) == 0:
                continue
            if header == None:
                header = line
                continue
            sequence_id = line[0]
            vquest_output[sequence_id] = {}
            for column,val in zip(header,line):
                vquest_output[sequence_id][column] = val
    return vquest_output

def get_assembly_alleles(assembly_allele_fasta):
    allele_names = []
    alleles = []
    assemble_alleles = list(SeqIO.parse(assembly_allele_fasta, "fasta"))
    assemble_alleles = [i for i in assemble_alleles if len(i.seq) >= 9]
    for i,assemble_allele in enumerate(assemble_alleles):
        allele_names.append(assemble_allele.id)
        assemble_allele.id = str(i)
        assemble_allele.name = ""
        assemble_allele.description = ""
        alleles.append(assemble_allele)
    return (allele_names,alleles)


parser = argparse.ArgumentParser(description='Process IGenotyper alleles')
parser.add_argument('assembly_genes',
                    help='assembly_genes.fasta IGenotyper output')
parser.add_argument('ccs_genes',
                    help='ccs_genes.fasta IGenotyper output')
parser.add_argument('samples_fofn',
                    help='A single column file with file paths of assemble_genes.fasta for other samples')
parser.add_argument('airr_vquest_output',
                    help='AIRR IMGT v-quest output file')
parser.add_argument('--species',  default="human",
                    help='Species. Must be an IMGT option')
parser.add_argument('--locus',  default="TR",
                    help='Receptor or Locus type. Must be an IMGT option')

args = parser.parse_args()

ccs_allele_fasta = ccs_read_count(args.ccs_genes)
other_samples_seqs = read_other_samples(args.samples_fofn)

allele_names,assemble_alleles = get_assembly_alleles(args.assembly_genes)

run_config = set_run_config(species=args.species,receptorOrLocusType=args.locus,fileSequences=args.assembly_genes)

run_vquest(run_config,assemble_alleles,args.airr_vquest_output)
vquest_output = parse_vquest(args.airr_vquest_output)

for allele_name,assemble_allele in zip(allele_names,assemble_alleles):    
    allele = ValidateAllele(allele_name,str(assemble_allele.seq))
    allele.add_imgt_data(vquest_output[assemble_allele.id])
    allele.add_ccs_support(ccs_allele_fasta)
    allele.found_in_other_samples(other_samples_seqs)
    print(allele)
    
