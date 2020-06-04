#!/usr/bin/env python3
############################################################
## Jose F. Sanchez                                        ##
## Copyright (C) 2020 Lauro Sumoy Lab, IGTP, Spain        ##
############################################################
"""
Calls IslandPath Dimob software for identification of putative genomics islands within a genome.
"""
## useful imports
import time
import io
import os
import re
import sys
import pandas as pd
import shutil
from termcolor import colored

## import my modules
from BacterialTyper.scripts import functions
from BacterialTyper.config import set_config

## to fix:  dimob_perl_script = "/imppc/labs/lslab/jsanchez/git_repo/BacterialTyper/third_party/IslandPath-DIMOB/Dimob.pl"

######
def GI_module(genbank_file, name, outdir, cutoff_dinuc_bias=8, min_length=1000, Debug):
    """Identify genomic islands (GI) within the genbank file provided. They are calculated
    based on gene annotation and dinucleotide bias region using the software `IslandPath-DIMOB`_.
    
    :param genbank_file: Absolute path to annotation file in Genbank format.
    :param name: Sample identifier. 
    :param outdir: Absolute path to output folder.
    :param cutoff_dinuc_bias: Dinucleotide bias cutoff
    :param min_length: Minimun length for the regions to be reported

    :type name: string
    :type genbank_file: string
    :type outdir: string
    :type cutoff_dinuc_bias: int
    :type min_length: int

    The Dimob.pl perl script has two mandatory argument which are the input :file:`genbank_file` and an output name.
    
    .. code-block:: sh

        Usage:
        perl Dimob.pl <genome.gbk> <output_name> [cutoff_dinuc_bias] [min_length]
        
        Default values:
            cutoff_dinuc_bias = 8
            min_length = 8000
        
        Example:
            perl Dimob.pl example/NC_003210.gbk NC_003210_GIs
            perl Dimob.pl example/NC_003210.gbk NC_003210_GIs 6 10000
            perl Dimob.pl example/NC_000913.embl NC_000913_GIs 6 10000

    
    During the development of BacterialTyper, we generated a modification of the original `IslandPath-DIMOB`_ to analyze 
    contig sequence data and generated different output format for better clarificaiton and interpretaion of results. 
    We forked the original code into a new git repository and update the code accordingly. See details here: https://github.com/JFsanchezherrero/islandpath.
    
     .. include:: ../../links.inc
    
    """
    
    ## filename stamp of the process
    filename_stamp = outdir + '/.Dimob'

    # check if previously done
    if os.path.isfile(filename_stamp):
        stamp =    functions.read_time_stamp(filename_stamp)
        print (colored("\tA previous command generated results on: %s [%s -- Dimob]" %(stamp, name), 'yellow'))
    else:    
        ## debug message
        if (Debug):
            print (colored("**DEBUG: Call Dimob for sample %s " %gbk_file + "**", 'yellow'))

        ## Call IslandPath Dimob executable perl file.
        dimob_pl = set_config.get_exe("dimob")
        perl_exe = set_config.get_exe("perl")
        perl_cmd = '%s %s %s %s %s %s' %(perl_exe, dimob_pl, genbank_file, outdir, cutoff_dinuc_bias, min_length)
    
        code = functions.system_call(perl_cmd)
        ##
        if code:
            ## when finished print time stamp in  output + '/.Dimob'
            stamp =    functions.print_time_stamp(filename_stamp)
        else:
            return False
    
    return (outdir)

def help_Dimob():
    print (colored("\n\n***** TODO: Generate this help message *****\n\n", 'red'))
    