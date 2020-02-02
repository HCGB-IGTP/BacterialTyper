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
from BacterialTyper import functions
from BacterialTyper import config

## to fix: 
dimob_perl_script = "/imppc/labs/lslab/jsanchez/git_repo/BacterialTyper/third_party/IslandPath-DIMOB/Dimob.pl"

def GI_module(genbank_file, name, outdir, cutoff_dinuc_bias=8, min_length=1000):
    """Identify genomic islands (GI) within the genbank file provided. They are calculated
    based on gene annotation and dinucleotide bias region using the software `IslandPath-DIMOB`_.
    
    :param genbank_file: Absolute path to annotation file in Genbank format.
    :type genbank_file: string
    
    :param name: Sample identifier. 
    :type name: string
    
    :param outdir: Absolute path to output folder.
    :type outdir: string
    
    :param cutoff_dinuc_bias: Dinucleotide bias cutoff
    :type cutoff_dinuc_bias: int
    
    :param min_length: Minimun length for the regions to be reported
    :type min_length: int
    
    
    The Dimob.pl perl script has two mandatory argument which are the input :file:`genbank_file` and an output name.
    
    .. code-block:: sh

        perl BacterialTyper/../third_party/IslandPath-DIMOB/Dimob.pl genbank_file output
    
    During the development of BacterialTyper, we generated a modification of the original `IslandPath-DIMOB`_ to analyze 
    contig sequence data and generated different output format for better clarificaiton and interpretaion of results. 
    We forked the original code into a new git repository and update the code accordingly. See details here: https://github.com/JFsanchezherrero/islandpath.
    
    .. seealso:: This function depends on other BacterialTyper functions called:
    
        - :func:`BacterialTyper.genomic_island.call_Dimob`
    
    
     .. include:: ../../links.inc
    
    """   
    
    

def call_Dimob(genbank_file, name, outdir, cutoff_dinuc_bias=8, min_length=8000):
    """Call IslandPath Dimob executable perl file.
    
    """