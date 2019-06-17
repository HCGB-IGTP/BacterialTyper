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
	
	### symbolic links
	print ("+ Retrieve all genomes assembled...")
	pd_samples_retrieved = sample_prepare.get_files(options, input_dir, "assembly")

	## annotate
	print ("+ Annotate assemblies using prokka:")
	print ("\t-Option: kingdom = ", options.kingdom,"; Annotation mode")
	if options.genera == 'Other':
		print ("\t-Option: genera = Off; No genus-specific BLAST databases option provided")
	else:
		print ("\t-Option: genera = ", options.genera,"; Genus-specific BLAST databases option provided")

	print ("\t-Option: addgenes; Add 'gene' features for each 'CDS' feature")
	print ("\t-Option: addmrna;  Add 'mRNA' features for each 'CDS' feature")
	print ("\t-Option: cdsrnaolap;  Allow [tr]RNA to overlap CDS")

	## optimize threads
	threads_module = functions.optimize_threads(options.threads, pd_samples_retrieved.index.size)

	## send for each sample
	with concurrent.futures.ThreadPoolExecutor(max_workers=int(options.threads)) as executor:
		commandsSent = { executor.submit(annotation.module_call, row['assembly'], options.kingdom, options.genera, outdir + '/' + row['samples'], row['samples'] , threads_module): index for index, row in pd_samples_retrieved.iterrows() }
		for cmd2 in concurrent.futures.as_completed(commandsSent):
			details = commandsSent[cmd2]
			try:
				data = cmd2.result()
			except Exception as exc:
				print ('***ERROR:')
				print (cmd2)
				print('%r generated an exception: %s' % (details, exc))

	## time stamp
	start_time_partial = functions.timestamp(start_time_total)
	
	exit()
	
	## Check each annotation using BUSCO
	functions.boxymcboxface("BUSCO Annotation Quality check")
	database_folder = os.path.abspath(options.database)
	outdir_BUSCO = functions.create_subfolder("BUSCO", outdir)
	BUSCO_Database = database_folder + '/BUSCO'
	dataFrame_results = qc.BUSCO_call(options.BUSCO_dbs, my_dir_annot, BUSCO_Database, outdir_BUSCO, options.threads, "protein")
	
	## functions.timestamp
	print ("+ Quality control of all samples finished: ")
	start_time_partial = functions.timestamp(start_time_partial_BUSCO)
	
	## summarize
	if (options.skip_report):
		print ("+ No report generation...")
	else:
		## generate plots
		print ("+ Generate summarizing plots...")
		qc.BUSCO_plots(dataFrame_results, outdir, options.threads)	
	
	## multiqc report plot
	## busco plot of samples or datasets per sample

	print ("\n*************** Finish *******************")
	start_time_partial = functions.timestamp(start_time_total)

	print ("+ Exiting Assembly module.")
	exit()

