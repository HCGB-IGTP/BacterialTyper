#!/usr/bin/env python3
##########################################################
## Jose F. Sanchez										##
## Copyright (C) 2019 Lauro Sumoy Lab, IGTP, Spain		##
##########################################################
'''
Calls card-trick module to parse CARD resistance information.
'''
## useful imports
import time
import io
import os
import re
import sys
import pandas as pd
from sys import argv
from io import open
from termcolor import colored
import card_trick

## import my modules
import HCGB
import HCGB.functions.time_functions as HCGB_time
import HCGB.functions.files_functions as HCGB_files
import HCGB.functions.main_functions as HCGB_main

from BacterialTyper.config import set_config

##########
def get_info_CARD(IDs, term, dataF):
	## dataF contains CARD ontology: downloaded and parse using prepare_card_data
	## IDs is an input list
	## term is the type of search to do using card_trick
	
	# search for terms provided	
	matching_terms = card_trick.ontology_functions.search(IDs, dataF, term, False)
	return (matching_terms)

#####################
def prepare_card_data(database_folder):
	
	## create CARD folder
	abs_folder = os.path.abspath(database_folder)
	CARD_folder = HCGB_files.create_subfolder('CARD', abs_folder)
	
	## make stamp time
	filename_stamp = CARD_folder + '/.success'

	if os.path.isfile(filename_stamp):
		stamp =	HCGB_time.read_time_stamp(filename_stamp)
		print (colored("\tA previous command generated results on: %s [CARD Ontology Data]" %stamp, 'yellow'))

		## check time passed
		days_passed = HCGB_time.get_diff_time(filename_stamp)
		print ("\t** %s days ago" %days_passed)		
		if (days_passed > 30): ## download again
			print ("\t ** Downloading information again just to be sure...")
			download=True
		else:
			print ("\t ** No need to download data again.")
			download=False
	else:
		download=True

	###
	if download:
		## uptade database in a path
		aro_obo_file = card_trick.ontology_functions.update_ontology(CARD_folder, False)
	
		## get ontology and save it in csv
		return_frame = card_trick.ontology_functions.parse_ontology(aro_obo_file, False)
	
		### if success return folder name
		if not return_frame.empty:
			## success stamps
			filename_stamp = CARD_folder + '/.success'
			stamp =	HCGB_time.print_time_stamp(filename_stamp)	
		else:
			return (FAIL)

	## return folder name
	return(CARD_folder)


def dead_code():
	## card_prepareref
	rename_info = card_prepareref + '00.rename_info'
	
	## outfile
	outfile = card_prepareref + '00.info_dictionary'
	out_file_handle = open(outfile, 'w')	
	
	## get info
	lines = HCGB_main.readList_fromFile(rename_info)
	for l in lines:
		names = l.split('\t')
		## original name \t ariba_name
		out_file_handle.write(names[1].split('.')[0] + '\t' + names[0].split('.')[0] + '\n')

	out_file_handle.close()	
	return (outfile)
	

