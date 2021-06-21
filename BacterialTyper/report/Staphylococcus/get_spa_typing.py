#!/usr/bin/env python3
##########################################################
## Jose F. Sanchez                                      ##
## Copyright (C) 2019-2020 Lauro Sumoy Lab, IGTP, Spain ##
##########################################################
from _ast import If
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
import pandas as pd

## import my modules
from BacterialTyper.config import set_config

## import my HCGB module 
import HCGB.functions.files_functions as HCGB_files
import HCGB.functions.aesthetics_functions as HCGB_aes
import HCGB.functions.main_functions as HCGB_main
import HCGB.functions.time_functions as HCGB_time

##############
def help_options():
    print ("\nUSAGE: python %s fasta database_folder...\n"  %os.path.realpath(__file__))

##############
def help_spaTyper():
    print (spaTyper.utils.extra_info())
    exit()
    
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
    
    timeStamp = HCGB_main.retrieve_matching_files(db_folder, 'success')
    
    ## info2return
    info_dict = { 'repeats' : spaTyper_repeats,
                 'url_repeats': "http://spa.ridom.de/dynamic/sparepeats.fasta",
                 'types': spaTyper_types,
                 'url_types': "http://spa.ridom.de/dynamic/spatypes.txt",
                 'folder': db_folder,
                 'stamp': HCGB_main.get_info_file(timeStamp[0])[0]
        }
    
    return(spaTyper_repeats, spaTyper_types, info_dict)

##############
def module_call(db_folder, dictionary_fasta_files, outdir_dict, debug):
    """
    
    """
        
    HCGB_files.create_folder(db_folder)
    if db_folder.endswith("spaTyper"):
        spaTyper_db = db_folder
    else:
        spaTyper_db = os.path.join(db_folder, "spaTyper")
        HCGB_files.create_folder(spaTyper_db)
        
    ## check if files are available
    (spaTyper_repeats, spaTyper_types, info_dict) = check_files(spaTyper_db, debug)

    ## Get the SpaTypes in fasta sequences
    seqDict, letDict, typeDict, seqLengths = spaTyper.spa_typing.getSpaTypes(spaTyper_repeats, spaTyper_types, debug)
    
    ## add info to return
    info_dict["letDict"] = letDict
    info_dict["seqLengths"] = list(seqLengths)
    
    ## debug messages
    if debug:
        HCGB_aes.debug_message("seqDict: Too large to print: See repeat_file for details", 'yellow')
        HCGB_aes.debug_message("typeDict: Too large to print: See repeat_order_file for details", 'yellow')
        HCGB_aes.debug_message("letDict: conversion dictionary", 'yellow')
        print (letDict)
        HCGB_aes.debug_message("seqLengths:", 'yellow')
        print (seqLengths)
        
    ## summary results
    results_summary = pd.DataFrame(columns=("sample", "sequence", "Repeats", "Repeat Type"))
    
    print("\n+ Create Staphylococcus Protein A (spa) typing for each sample:")
    
    ## for each sample get spaType
    for key, value in dictionary_fasta_files.items():
        
        ## save in folder for each sample
        report_folder = HCGB_files.create_folder(outdir_dict[key])
        spaType_folder = HCGB_files.create_subfolder('spatype', report_folder)
        results_file = os.path.join(spaType_folder, 'spatyper_results.txt')    
        stampfile = os.path.join(spaType_folder, '.success')
        ## check if previously done
        if os.path.isfile(stampfile):
            stamp = HCGB_time.read_time_stamp(stampfile)
            print (colored("\tA previous command generated results on: %s [%s]" %(stamp, key), 'yellow'))
            results_spa = HCGB_main.get_info_file(results_file)
            for i in results_spa:
                i_list = i.split(';')
                
                results_summary.loc[len(results_summary)] = (key, i_list[0].split(":")[1], 
                                                             i_list[1].split(":")[1], 
                                                             i_list[2].split(":")[1])    
        else:
            print ("\t+ Sample: ", key)
            returned_value = call_spaTyper(value, seqDict, letDict, typeDict, seqLengths, debug)
            
            if len(returned_value.keys()) > 1:
                print (colored("** Attention: >1 spaTypes detected for sample: %s" %key, 'red'))

            list_results = []
            for j in returned_value.keys():
                splitted = returned_value[j].split('::')
                results_summary.loc[len(results_summary)] = (key, j, splitted[2], splitted[1])    
                res_string = "Sequence name: " + j +  "; Repeats: " + splitted[2] + "; Repeat Type: " +  splitted[1]
                
                ## save into file
                list_results.append(res_string)
                
                ## debug messages
                if debug:
                    HCGB_aes.debug_message(res_string, "yellow")
                
            ## dump results in file
            HCGB_main.printList2file(results_file, list_results)
            
            ## print time stamp
            ## dump results in file
            HCGB_time.print_time_stamp(stampfile)
            
    ##
    return (results_summary, info_dict)

##############
def call_spaTyper(fasta_file, seqDict, letDict, typeDict, seqLengths, debug):
    """
    Call spaTyper for a fasta file provided using precomputed spa repeats orders and types.
    
    :param fasta_file: Assembly fasta file to check for spa repeats.
    :param seqDict: Sequence dictionary
    :param letDict:
    :param typeDict:
    :param seqLengths:
    :param debug:
    
    :type fasta_file:
    :type seqDict:
    :type letDict:
    :type typeDict:
    :type seqLengths:
    :type debug:
    """
    
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
    debug=False

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