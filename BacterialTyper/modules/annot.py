#!/usr/bin/env python3
'''
This code calls Prokka to generate a functional assembly annotation
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
from BacterialTyper import BUSCO_caller
from BacterialTyper.modules import qc

from BacterialTyper.modules import sample_prepare

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
	functions.boxymcboxface("Assembly annotation")

	print ("--------- Starting Process ---------")
	functions.print_time()

	## print further information for Prokka	
	if (options.help_Prokka):
		annotation.print_list_prokka()
		exit()
	
	## print further information for BUSCO databases	
	if (options.help_BUSCO):
		BUSCO_caller.print_help_BUSCO()
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
	
	## annotate
	my_dir_annot = {}
	for assemblies in my_assembly_list:
		sample_dir = outdir + '/sample'
		name = 'sample'
		dir_annot = annotation.module_call(sequence_fasta, options.kingdom, sample_dir, name, options.threads)	
		my_dir_annot[name] = dir_annot
	
	## Check each annotation using BUSCO
	functions.boxymcboxface("BUSCO Annotation Quality check")
	database_folder = os.path.abspath(options.database)
	outdir_BUSCO = functions.create_subfolder("BUSCO", outdir)
	BUSCO_Database = database_folder + '/BUSCO'
	dataFrame_results = qc.BUSCO_call(options.BUSCO_dbs, my_dir_annot, BUSCO_Database, outdir_BUSCO, options.threads, "protein")
	
	## functions.timestamp
	print ("+ Quality control of all samples finished: ")
	start_time_partial = functions.timestamp(start_time_partial_BUSCO)
	
	## generate plots
	print ("+ Generate summarizing plots...")
	qc.BUSCO_plots(dataFrame_results, outdir, options.threads)	
	
	## multiqc report plot
	## busco plot of samples or datasets per sample

	print ("\n*************** Finish *******************")
	start_time_partial = functions.timestamp(start_time_total)

	print ("+ Exiting Assembly module.")
	exit()
		
