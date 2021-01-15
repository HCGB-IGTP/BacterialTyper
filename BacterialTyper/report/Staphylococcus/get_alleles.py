#!/usr/bin/env python3
##########################################################
## Jose F. Sanchez                                      ##
## Copyright (C) 2019-2020 Lauro Sumoy Lab, IGTP, Spain ##
##########################################################
from builtins import print
"""
Searches for alleles sequences provided
"""
## useful imports
import time
import io
import os
import re
import sys
from sys import argv
from io import open
from termcolor import colored
import pandas as pd
from Bio import SeqIO

## import my modules
#from BacterialTyper.scripts import functions
from BacterialTyper.config import set_config
from BacterialTyper.scripts import edirect_caller

## import my HCGB module 
from HCGB.functions import files_functions
from HCGB.functions import time_functions
from HCGB.functions import main_functions

##############
def help_options():
    print ("\nUSAGE: python %s csv_file_alleles path sample...\n"  %os.path.realpath(__file__))

##############
def get_sequence(path, probe_ID, nucc_entry, format, revcomp=False, start=0, stop=-1):
    ## nucc_entry == e.g. NZ_CP029680.1 Staphylococcus aureus strain AR_0215 chromosome, complete genome
    
    ## create name
    nucc_entry_name = ""
    if (start != 0):
        nucc_entry_name = nucc_entry + "-start_" + str(start) + "-end_" + str(stop)
    else:
        nucc_entry_name = nucc_entry + "-all"
    
    if revcomp:
        nucc_entry_name = nucc_entry_name + "_revcomp"
    
    ## fasta file
    out_fasta_file = os.path.join(path, nucc_entry_name + ".fasta")

    ## check if previously retrieved
    filename_stamp = path + '/.success_' + nucc_entry_name
    if os.path.isfile(filename_stamp) and files_functions.is_non_zero_file(out_fasta_file):
        stamp = time_functions.read_time_stamp(filename_stamp)
        print (colored("\tA previous command generated results on: %s [%s -- %s]" %(stamp, probe_ID, nucc_entry_name), 'yellow'))
    else: 
        edirect_caller.generate_seq_search_call('nuccore', nucc_entry, out_fasta_file, revcomp, start, stop, format)
        ## timestamp
        stamp = time_functions.print_time_stamp(filename_stamp)
    
    for record in SeqIO.parse(out_fasta_file, 'fasta'):
        return(str(record.seq), record.id)


def get_alleles(csv_file, path, debug):
    
    ## get probe information
    probe_df = main_functions.get_data(csv_file, ',', 'index_col=0')
    
    ## debugging messages
    if debug:
        print ("### DEBUG: probe_df ###")
        print (probe_df)
        print ("### DEBUG: probe_df ###")
        
    ## This dataframe can contain several columns some of which are optional:
    ## - ID: Probe ID [mandatory]
    ## - Gene: Gene name [mandatory]
    ## - Variant: Variant name [mandatory]
    ## - Probe_Name: Probe name [optional]

    ## Probe information [mandatory]:
    ## Provide either the NCBI id and coordinates or the sequence of the probe to search
    ## - Probe_Coord: Provide the NCBI id and coordinates within square brackets. e.g. AF288215.1[1255:1275]
    ## - Probe_Seq: Probe sequence to search
    
    ## Primer information [optional]
    ## Provide if desired a primer to double check. Provide either the NCBI id and coordinates 
    ## or the primer sequence to search
    ## - Primer_Coord: Same format as for Probe_Coord, e.g. AF288215.1[1255:1275:r]
    ## - Primer_Seq: Sequence to search
    
    ## groupby variant
    variant_frame = probe_df.groupby(["Variant"])
    for name, frame in variant_frame:
        ## debugging messages
        if debug:
            print ("### DEBUG: frame ###")
            print (frame)
            print()

        ## search each probe
        for index, row in frame.iterrows():
            if (row['Probe_Coord']):
                (nucc_entry, start, stop, revcomp) = parse_probe(row['Probe_Coord'])
                ## debugging messages
                if debug:
                    print ("### DEBUG: parse_probe results ###")
                    print ("ID: " + index)
                    print ("Probe_Coord")
                    print ("nucc_entry: " + nucc_entry)
                    print ("start: " + start)
                    print ("stop: " + stop)
                    print ("revcomp: " + str(revcomp))
                
                ## sequence is save in file
                (probe_fasta, probe_fasta_header) = get_sequence(path, index+ '::Probe', nucc_entry, "fasta", revcomp, start, stop)
            else:
                if (row['Probe_Seq']):
                    probe_fasta = row['Probe_Seq']
                else:
                    print ("** ERROR: This probe does not contain either coordinates or sequence information.")
                    print (row)
                    exit()
            
            if (row['Primer_Coord']):
                (nucc_entry, start, stop, revcomp) = parse_probe(row['Primer_Coord'])
                ## debugging messages
                if debug:
                    print ("### DEBUG: parse_probe results ###")
                    print ("ID: " + index)
                    print ("Primer_Coord")
                    print ("nucc_entry: " + nucc_entry)
                    print ("start: " + start)
                    print ("stop: " + stop)
                    print ("revcomp: " + str(revcomp))
                
                ## sequence is save in file
                (primer_fasta, primer_fasta_header) = get_sequence(path, index + '::Primer', nucc_entry, "fasta", revcomp, start, stop)
            else:
                if (row['Primer_Seq']):
                    primer_fasta = row['Primer_Seq']
        
            ## generate results for each
            
            
            
                
def parse_probe(string_probe):
    ##
    string_search = re.search(r"(.*)\[(.*)\]", string_probe)
    split_brackets = string_search.group(2).split(":")
    
    ## e.g. AF026120.1[3:29]
    ## AF026120.1 # string_search.group(1)
    ## ['3', '29'] split_brackets
    
    ## reverse complement
    revComp = False
    if len(split_brackets) > 2:
        if split_brackets[2] == 'r':
            revComp = True
    
    ## return
    return (string_search.group(1), split_brackets[0], split_brackets[1], revComp)
    

def main():

    ## control if options provided or help
    if len(sys.argv) > 1:
        print ("")
    else:
        help_options()
        exit()
    
    ##
    csv_file = os.path.abspath(sys.argv[1])
    path = os.path.abspath(sys.argv[2])

    ##
    get_alleles(csv_file, path, False)
    

'''******************************************'''
if __name__== "__main__":
    main()

