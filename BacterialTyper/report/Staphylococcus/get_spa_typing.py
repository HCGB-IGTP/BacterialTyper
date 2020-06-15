#!/usr/bin/env python3
##########################################################
## Jose F. Sanchez                                      ##
## Copyright (C) 2019-2020 Lauro Sumoy Lab, IGTP, Spain ##
##########################################################
"""
Calls spaTyper for obtaining Staphylococcus aureus Spa protein type.
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
import spaTyper

## import my modules
from BacterialTyper.scripts import functions
from BacterialTyper.config import set_config


##############
def help_options():
    print ("\nUSAGE: python %s fasta database_folder...\n"  %os.path.realpath(__file__))

##############
def help_spaTyper():
    print (colored("\n\n***** TODO: Generate this help message *****\n\n", 'red'))

##############
def check_files(db_folder, debug):
    """
    Check if sparepeats and spatypes files are available.
    
    If not available, it downloads them from SeqNet/Ridom Spa Server
    
    :param db_folder: Database containing files
    :param debug: True/false for debugging messages
    
    :type db_folder: string
    :type debug: boolean 
    """
    
    spaTyper_repeats = os.path.join(db_folder, "sparepeats.fasta")
    spaTyper_types = os.path.join(db_folder, "spatypes.txt")
    
    print ("+ Check or downloads repeats fasta and types file in folder: ", db_folder)
    
    print ("Check file: http://spa.ridom.de/dynamic/sparepeats.fasta")
    spaTyper_repeats = spaTyper.utils.download_file_repeats(db_folder, debug)
    
    print ("Check file: http://spa.ridom.de/dynamic/spatypes.txt")
    spaTyper_types = spaTyper.utils.download_file_types(db_folder, debug)
        
    return(spaTyper_repeats, spaTyper_types)

##############
def module_call(db_folder, dictionary_fasta_files, debug):
    
    functions.create_folder(db_folder)
    if db_folder.endswith("spaTyper"):
        spaTyper_db = db_folder
    else:
        spaTyper_db = os.path.join(db_folder, "spaTyper")
        functions.create_folder(spaTyper_db)
        
    ## check if files are available
    (spaTyper_repeats, spaTyper_types) = check_files(spaTyper_db, debug)

    ## Get the SpaTypes in fasta sequences
    seqDict, letDict, typeDict, seqLengths = spaTyper.spa_typing.getSpaTypes(spaTyper_repeats, spaTyper_types, debug)
    
    ## debug messages
    if args.debug:
        print ('## Debug: seqDict: Too large to print: See repeat_file for details')
        print ('## Debug: typeDict: Too large to print: See repeat_order_file for details')
        print ('## Debug: letDict: conversion dictionary')
        print (letDict)
        print ('## Debug: seqLengths:')
        print (seqLengths)
        
    ## summary results
    results_summary = pd.DataFrame(columns=("sample", "sequence", "Repeats", "Repeat Type"))
    
    ## for each sample get spaType
    for key, value in dictionary_fasta_files.items():
        print ("+ Sample: ", key)
        returned_value = call_spaTyper(value, qDict, seqDict, letDict, typeDict, seqLengths, debug)
        
        for j in returned_value.keys():
            splitted = returned_value[j].split('::')
            if len(returned_value.keys()) > 1:
                print (colored("** Attention: >1 spaTypes detected for sample: %s" %key, 'red'))
            
            results_summary.loc[len(results_summary)] = (key, j, splitted[2], splitted[1])    
            print("Sequence name: ",j, "Repeats:", splitted[2], "Repeat Type:", splitted[1], '\n')    
    
    ##
    return (results_summary)

##############
def call_spaTyper(fasta_file, seqDict, letDict, typeDict, seqLengths, debug):
    
    #######################
    ## findPatterns for each fasta file
    #######################
    ## get sequence dictionary
    qDict = spaTyper.utils.fasta_dict(fasta_file)
    
    ## find pattern
    do_enrich = False
    the_out = spaTyper.spa_typing.findPattern(qDict, seqDict, letDict, typeDict, seqLengths, do_enrich, debug)
    return the_out

##############
def main():

    ## control if options provided or help
    if len(sys.argv) > 1:
        print ("")
    else:
        help_options()
        exit()
        
    fasta_file = os.path.abspath(argv[1])
    db_folder = os.path.abspath(argv[2])
    debug=True

    ## get pattern
    (spaTyper_repeats, spaTyper_types) = check_files(db_folder, debug)
    
    ## Get the SpaTypes in fasta sequences
    seqDict, letDict, typeDict, seqLengths = spaTyper.spa_typing.getSpaTypes(spaTyper_repeats, spaTyper_types, debug)
    
    ## debug messages
    if debug:
        print ('## Debug: seqDict: Too large to print: See repeat_file for details')
        print ()
        #print (seqDict)
        print ('## Debug: letDict: conversion dictionary')
        print (letDict)
        print ()
        print ('## Debug: typeDict: Too large to print: See repeat_order_file for details')
        #print (typeDict)
        print ()
        print ('## Debug: seqLengths:')
        print (seqLengths)
        print ()
    
    returned_value = call_spaTyper(fasta_file, seqDict, letDict, typeDict, seqLengths, debug)
    
    if len(returned_value.keys()) > 1:
        print ("** Attention: >1 spaTypes detected")
    
    for j in returned_value.keys():
        splitted = returned_value[j].split('::')
        print("Sequence name: ",j, "Repeats:", splitted[2], "Repeat Type:", splitted[1], '\n')    

    

'''******************************************'''
if __name__== "__main__":
    main()