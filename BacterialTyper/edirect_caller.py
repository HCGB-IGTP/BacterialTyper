#usr/bin/en python3
'''
This code...
Jose F. Sanchez
Copyright (C) 2019 Lauro Sumoy Lab"," IGTP"," Spain
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
from BacterialTyper import functions
from BacterialTyper import config

####################################
def help_edirect():
	functions.boxymcboxface("ENTREZ DIRECT (EDirect)")
	print ("+ Connects with NCBI entrez databases")
	print ("+ Searches, retrieves, and parses data from NCBI databases through the Unix command line")
	print ("+ Provides access to the NCBI's suite of interconnected databases")
	print ("+ Includes an argument-driven function that simplifies the extraction of \ndata from document summaries or other results that are in structured XML format.")
	print ("\n")

###############
def generate_docsum_call(db, query, outfile):
	esearch_bin = config.get_exe("esearch") 
	efetch_bin = config.get_exe("efetch")
	return(docsum_call(db, query, outfile, esearch_bin, efetch_bin))

###############
def docsum_call(db, query, outfile, esearch_bin, efetch_bin):
	cmd = ("%s -db %s -query %s | %s -format docsum > %s" %(esearch_bin, db, query, efetch_bin, outfile))
	return(functions.system_call(cmd))
	
###############
def generate_xtract_call(docsum_file, pattern, element, outfile):
	xtract_bin = config.get_exe("xtract") 
	return(xtract_call(docsum_file, pattern, element, outfile, xtract_bin))
	
###############
def xtract_call(docsum_file, pattern, element, outfile, xtract_bin):
	cmd = ("cat %s | %s -pattern %s -element %s > %s" %(docsum_file, xtract_bin, pattern, element, outfile))
	return(functions.system_call(cmd))	

###############
def main():
	## this code runs when call as a single script

  	## control if options provided or help
	if len(sys.argv) > 1:
		print ("")
	else:
		help_edirect()
		exit()

	### TODO implement main function
	## to call independently

######
'''******************************************'''
if __name__== "__main__":
	main()
		
	

## docsum Nucleotide entry
##esearch -db nuccore -query NZ_CP029083.1 | efetch -format docsum > docsum_nucleotide_entry.txt
##cat docsum_nucleotide_entry.txt | xtract -pattern DocumentSummary -element BioSample > BioSample.txt

## Docsum BioSample
##esearch -query BioSample.txt -db assembly | efetch -format docsum > BioSample_docsum.txt
##cat BioSample_docsum.txt | xtract -pattern DocumentSummary -element FtpPath_GenBank
