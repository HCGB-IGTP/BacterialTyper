#!/usr/bin/env python3
############################################################
## Jose F. Sanchez                                        ##
## Copyright (C) 2019-2021 Lauro Sumoy Lab, IGTP, Spain   ##
############################################################
from pandas.tests.io.json.conftest import orient
"""
Calls assembly_stats pip module
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
import assembly_stats

##
import HCGB.functions.aesthetics_functions as HCGB_aes
import HCGB.functions.main_functions as HCGB_main

############
def assembly_stats_caller(fasta_file, out_file, debug):
    
    contig_lens, scaffold_lens, gc_cont = assembly_stats.read_genome(fasta_file)
    
    ## debug messages
    if debug:
        HCGB_aes.debug_message("contig_lens", "yellow")
        print(contig_lens)
        HCGB_aes.debug_message("scaffold_lens", "yellow")
        print(scaffold_lens)
        HCGB_aes.debug_message("gc_cont", "yellow")
        print(gc_cont)
    
    ## get stats
    contig_stats = assembly_stats.calculate_stats(contig_lens, gc_cont)
    scaffold_stats = assembly_stats.calculate_stats(scaffold_lens, gc_cont)
    
    ## debug messages
    if debug:
        HCGB_aes.debug_message("contig_stats", "yellow")
        print(contig_stats)
        HCGB_aes.debug_message("scaffold_stats", "yellow")
        print(scaffold_stats)
    
    stat_output = {'Contig Stats': contig_stats,
                   'Scaffold Stats': scaffold_stats}
    
    ## save results in file
    HCGB_main.printDict2file(out_file + '-contigs.csv', contig_stats, ",")
    HCGB_main.printDict2file(out_file + '-scaffolds.csv', scaffold_stats, ",")
    
    ## create stats in excel file
    assembly_stats_file = out_file + '_stats.xlsx'
    parse_stats(stat_output, assembly_stats_file, debug)
    
    return(stat_output, assembly_stats_file)

####################################
def parse_stats(stat_output, assembly_stats_file, debug):
    
    ## create single excel file
    ## write to excel
    writer_Excel = pd.ExcelWriter(assembly_stats_file, engine='xlsxwriter') ## open excel handle
    
    ## loop over dictionary and print dataframe
    for name, each_data in stat_output.items():
        if debug:
            print ("Name: ", name)
            print ("Data: ", each_data)
            
        each_df = pd.DataFrame.from_dict(each_data, orient='index')
        each_df.to_excel(writer_Excel, sheet_name=name) ## write excel handle
            
    writer_Excel.save() ## close excel handle

############
def help_options():
    print ("\nUSAGE:\npython %s fasta_file out_file_name\n"  %os.path.realpath(__file__))

############
def main():
    ## this code runs when call as a single script

    ## control if options provided or help
    if len(sys.argv) > 1:
        print ("")
    else:
        help_options()
        exit()        

    file_in = os.path.abspath(argv[1])
    file_out = os.path.abspath(argv[2])
    
    ##
    path_to_sample = assembly_stats_caller(file_in, file_out, True)
        
############
'''******************************************'''
if __name__== "__main__":
    main()


