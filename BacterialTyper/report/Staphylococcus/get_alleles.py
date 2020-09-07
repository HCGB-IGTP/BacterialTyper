#!/usr/bin/env python3
##########################################################
## Jose F. Sanchez                                      ##
## Copyright (C) 2019-2020 Lauro Sumoy Lab, IGTP, Spain ##
##########################################################
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

## import my modules
from BacterialTyper.scripts import functions
from BacterialTyper.config import set_config
from BacterialTyper.scripts import edirect_caller

##############
def help_options():
    print ("\nUSAGE: python %s csv_file_alleles sample...\n"  %os.path.realpath(__file__))

##############
def get_sequence(path, nucc_entry, tmp_folder, format, revcomp=False, start=0, stop=-1):
    ## nucc_entry == e.g. NZ_CP029680.1 Staphylococcus aureus strain AR_0215 chromosome, complete genome
    
    ## create name
    nucc_entry_name = ""
    if (start != 0):
        print()
        nucc_entry_name = nucc_entry + "-start_" + str(start) + "-end_" + str(stop)
    else:
        nucc_entry_name = nucc_entry + "-all"
    
    if revcomp:
        nucc_entry_name = nucc_entry_name + "_revcomp"
    
    ##
    out_fasta_file = os.path.join(path, nucc_entry_name + ".fasta")
    edirect_caller.generate_seq_search_call('nuccore', nucc_entry, out_fasta_file, revcomp, start, stop, format)


##############
def main():

    ## control if options provided or help
    if len(sys.argv) > 1:
        print ("")
    else:
        help_options()
        exit()
        
        
    nucc_entry = sys.argv[1]
    tmp_folder = os.path.abspath(sys.argv[2])
    ## e.g. AB043554.1 [824:852]
    get_sequence(tmp_folder, nucc_entry, tmp_folder, 'fasta', revcomp=True, start=sys.argv[3], stop=sys.argv[4])
    

'''******************************************'''
if __name__== "__main__":
    main()

