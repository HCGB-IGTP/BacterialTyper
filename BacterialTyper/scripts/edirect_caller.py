#!/usr/bin/env python3
##########################################################
## Jose F. Sanchez										##
## Copyright (C) 2019 Lauro Sumoy Lab, IGTP, Spain		##
##########################################################
'''
Calls Entrez Direct NCBI search.
'''
## useful imports
import time
import io
import os
import re
import sys
import pandas as pd
from termcolor import colored

## import my modules
import HCGB
import HCGB.functions.time_functions as HCGB_time
import HCGB.functions.system_call_functions as HCGB_sys
import HCGB.functions.aesthetics_functions as HCGB_aes
from BacterialTyper.config import set_config

## additional information:
## https://ncbi-hackathons.github.io/EDirectCookbook/

####################################
def help_edirect():
	"""
	Prints edirect_caller script help options
	"""
	HCGB_aes.boxymcboxface("ENTREZ DIRECT (EDirect)")
	print ("+ Connects with NCBI entrez databases")
	print ("+ Searches, retrieves, and parses data from NCBI databases through the Unix command line")
	print ("+ Provides access to the NCBI's suite of interconnected databases")
	print ("+ Includes an argument-driven function that simplifies the extraction of \ndata from document summaries or other results that are in structured XML format.")
	print ("\n")

###############
def generate_docsum_call(db, query, outfile):
	esearch_bin = set_config.get_exe("esearch") 
	efetch_bin = set_config.get_exe("efetch")
	return(docsum_call(db, query, outfile, esearch_bin, efetch_bin))

###############
def docsum_call(db, query, outfile, esearch_bin, efetch_bin):
	cmd = ("%s -db %s -query %s | %s -format docsum > %s" %(esearch_bin, db, query, efetch_bin, outfile))
	return(HCGB_sys.system_call(cmd))
	
###############
def generate_xtract_call(docsum_file, pattern, element, outfile):
	xtract_bin = set_config.get_exe("xtract") 
	return(xtract_call(docsum_file, pattern, element, outfile, xtract_bin))
	
###############
def xtract_call(docsum_file, pattern, element, outfile, xtract_bin):
	cmd = ("cat %s | %s -pattern %s -sep ',' -element %s > %s" %(docsum_file, xtract_bin, pattern, element, outfile))
	return(HCGB_sys.system_call(cmd))	

###############
def generate_seq_search_call(db, query, outfile, revcomp, start=0, end=-1, format='fasta'):

	## Sequence Range
	##   -seq_start     First sequence position to retrieve
	##   -seq_stop      Last sequence position to retrieve
	##   -strand        1 = forward DNA strand, 2 = reverse complement
	##   -revcomp       Shortcut for strand 2

	efetch_bin = set_config.get_exe("efetch")
	cmd = ("%s -db %s -id %s -seq_start %s -seq_stop %s -format %s" %(efetch_bin, db, query, start, end, format))
	
	## add reverse complement
	if (revcomp):
		cmd = cmd + ' -revcomp'

	## add output file
	cmd = cmd + ' > %s' %outfile
		
	return(HCGB_sys.system_call(cmd))

## docsum Nucleotide entry
##esearch -db nuccore -query NZ_CP029083.1 | efetch -format docsum > docsum_nucleotide_entry.txt
##cat docsum_nucleotide_entry.txt | xtract -pattern DocumentSummary -element BioSample > BioSample.txt

## Docsum BioSample
##esearch -query BioSample.txt -db assembly | efetch -format docsum > BioSample_docsum.txt
##cat BioSample_docsum.txt | xtract -pattern DocumentSummary -element FtpPath_GenBank

