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
def retrieve_genes_ids_sequences(profile, gene_ID, debug):
    """    
    """
    ## given a profile folder
    if debug:
        print ("** DEBUG **")
        print ("profile: ", profile)
        print ("gene_id: ", gene_ID)
        
    ##
    assembled_genes_list = functions.main_functions.retrieve_matching_files(profile, "assembled_genes.fa", debug)
    assembled_genes_list = [s for s in assembled_genes_list if 'ariba.tmp' not in s]
    
    if debug:
        print ("** DEBUG **")
        print ("assembled_genes_list: ", assembled_genes_list)

    if os.path.isfile(assembled_genes_list[0]):
        for record in SeqIO.parse(assembled_genes_list[0], "fasta"):
            if debug:
                print ("** DEBUG: ", record.description, " **")
                #print (record.seq)
 
            search_ID = re.search(gene_ID, record.description)
            if (search_ID):
                return (record.id, str(record.seq))

        return('','')

## retrieve info
def get_genes_profile(samples_info, gene_names, debug):
    """    
    """
    results_profileIDs = pd.DataFrame(columns=('sample', 'gene', 'value'))
    sample_frame = samples_info.groupby(["name"])
    for g in gene_names:
        print ("\t+", g)
        for name, cluster_df in sample_frame:
            my_list_profiles = cluster_df.loc[cluster_df['tag'] == 'profile']['ext'].to_list()
            if debug:
                print ("name: ", name)
                print ("my_list_profiles:")
                print (my_list_profiles)
            
            for p in my_list_profiles:
                profile_csv = cluster_df.loc[cluster_df['ext'] == p]['sample'][0]
                value = retrieve_genes_ids_profile(profile_csv, g, debug)
                ## save results 
                if (value == 'no'):
                    results_profileIDs.loc[len(results_profileIDs)] = (name, g, value)
                else:
                    for colName, colData in value.iteritems():
                        results_profileIDs.loc[len(results_profileIDs)] = (name, colName, colData) 
    
    return (results_profileIDs)

##############
def retrieve_genes_ids_profile(profile, gene_ID, debug):
    """    
    """
    ## read data    
    get_csv_data = functions.main_functions.get_data(profile, ',', '')
    
    ## debug messages
    if debug:
        print ("** DEBUG **")
        print ("profile: ", profile)
        print ("gene_id: ", gene_ID)
        print ("data")
        print (get_csv_data)
    
    if (get_csv_data.filter(regex=('^' + gene_ID + '.*')).empty):
        return('no')
    else:
        this_df = get_csv_data.filter(regex=('^' + gene_ID + '.*'))
        filtered_genes = this_df.columns.values
        ## debug messages
        if debug:
            print ("** DEBUG **")
            print (filtered_genes)
            print (this_df[filtered_genes])
        
        return (this_df[filtered_genes])     
    
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
    
    retrieve_genes_ids_profile(profile, gene_ID, True)

'''******************************************'''
if __name__== "__main__":
    main()
