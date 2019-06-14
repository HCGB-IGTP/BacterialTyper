#!/usr/bin/env python3
'''
This code calls several programs to generate a Mobile Genetic Element analysis
Jose F. Sanchez
Copyright (C) 2019 Lauro Sumoy Lab, IGTP, Spain
'''
## useful imports
import time
import io
import os
import re
import sys
from io import open
import concurrent.futures
from termcolor import colored
import pandas as pd

## import my modules
from BacterialTyper import functions
from BacterialTyper import config
from BacterialTyper import annotation
from BacterialTyper import bacteriophage
from BacterialTyper.modules import sample_prepare
import PhiSpy

####################################
def run(options):

	## init time
	start_time_total = time.time()

	## debugging messages
	global Debug
	if (options.debug):
		Debug = True
	else:
		Debug = False

	### species_identification_KMA -> most similar taxa
	functions.pipeline_header()
	functions.boxymcboxface("Mobile Genetic Elements (MGE) identification")

	print ("--------- Starting Process ---------")
	functions.print_time()
	
	## print information and exit
	if options.help_PhiSpy:
		PhiSpy.print_list()
		exit()
		
	## absolute path for in & out
	input_dir = os.path.abspath(options.input)
	outdir = os.path.abspath(options.output_folder)

	## generate output folder
	functions.create_folder(outdir)

	### time stamp
	start_time_partial = start_time_total
	start_time_partial_assembly = start_time_partial
	
	### symbolic links
	print ("+ Retrieve all genomes assembled...")
	pd_samples_retrieved = sample_prepare.get_files(options, input_dir, "assembly")
	
	print (pd_samples_retrieved)
	exit()

	

