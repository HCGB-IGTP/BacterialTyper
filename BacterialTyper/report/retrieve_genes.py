#!/usr/bin/env python3
##########################################################
## Jose F. Sanchez                                      ##
## Copyright (C) 2019-2020 Lauro Sumoy Lab, IGTP, Spain ##
##########################################################
"""
Retrieves genes ids from profile analysis generated
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
from HCGB import functions
from BacterialTyper.config import set_config

##############
def help_options():
    print ("\nUSAGE: python %s gene_id profile_folder...\n"  %os.path.realpath(__file__))

##############
def help_retrieve_genes():
    print (colored("\n\n***** TODO: Generate this help message *****\n\n", 'red'))

##############
def retrieve_genes(profile, gene_ID, debug):
    ## given a profile folder
    if debug:
        print ("** DEBUG **")
        print ("profile: ", profile)
        print ("gene_id: ", gene_ID)
        
    ##
    assembled_genes = functions.main_functions.retrieve_matching_files(profile, "assembled_genes.fa", debug)[0]
    if os.path.isfile(assembled_genes):
        for record in SeqIO.parse(assembled_seqs, "fasta"):
            if debug:
                print ("** DEBUG **")
                print (record.description)
                print (record.seq)
            
            search_ID = re.search(gene_ID, record.description)
            if (search_ID):
                print (record.description)
                print (record.seq)
       
                return (record.description, record.seq)

        ## no return
        return("n.a", "n.a")        

##############
def main():

    ## control if options provided or help
    if len(sys.argv) > 1:
        print ("")
    else:
        help_options()
        exit()

    ##
    gene_ID = sys.argv[1]
    profile = os.path.abspath( sys.argv[2] )
    
    retrieve_genes(profile, gene_ID, False)

'''******************************************'''
if __name__== "__main__":
    main()