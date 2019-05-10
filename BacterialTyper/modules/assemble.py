#!/usr/bin/env python3
'''
This code calls different submodules to assemble each sample using SPADES, check quality using BUSCO and annotate using PROKKA
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

## import my modules
from BacterialTyper import functions
from BacterialTyper import config

from BacterialTyper import spades_assembler
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

	## message header
	functions.pipeline_header()
	functions.boxymcboxface("Assembly module")
	print ("--------- Starting Process ---------")
	functions.print_time()
	
	## print further information for BUSCO databases	
	if (options.help_BUSCO):
		BUSCO_caller.print_help_BUSCO()
		exit()

	## absolute path for in & out
	input_dir = os.path.abspath(options.input)
	outdir = os.path.abspath(options.output_folder)

	## get files
	pd_samples_retrieved = sample_prepare.get_files(options, input_dir)
	## generate output folder
	functions.create_folder(outdir)

	### call assemble using spades
	start_time_partial = start_time_total
	start_time_partial_assembly = start_time_partial
	
	## optimize threads
	workers = ""
	if (options.threads > 12): ## if many cpus provided...
		threads_module = functions.optimize_threads(options.threads, pd_samples_retrieved.index.size)
		workers = pd_samples_retrieved.index.size
	elif (options.threads > 6):
		threads_module = int(options.threads/2)
		workers = 2
	else:
		workers = 1
		threads_module = options.threads

	print ('+ Running modules SPADES...')

	# We can use a with statement to ensure threads are cleaned up promptly
	with concurrent.futures.ThreadPoolExecutor(max_workers=int(workers)) as executor:
		## send for each sample
		commandsSent = { executor.submit( check_sample_assembly, row['samples'], outdir + '/' + row['samples'],  row['R1'], row['R2'], threads_module): index for index, row in pd_samples_retrieved.iterrows() }

		for cmd2 in concurrent.futures.as_completed(commandsSent):
			details = commandsSent[cmd2]
			try:
				data = cmd2.result()
			except Exception as exc:
				print ('***ERROR:')
				print (cmd2)
				print('%r generated an exception: %s' % (details, exc))
		
	## functions.timestamp
	print ("Assembly of all samples finished: ")
	start_time_partial = functions.timestamp(start_time_partial_assembly)
	start_time_partial_BUSCO = start_time_partial

	## Check each 
	qc.assembly_check(options)
	
	## functions.timestamp
	print ("Quality control of all samples finished: ")
	start_time_partial = functions.timestamp(start_time_partial_BUSCO)

	print ("\n*************** Finish *******************")
	start_time_partial = functions.timestamp(start_time_total)

	print ("+ Exiting Assembly module.")
	exit()

def check_sample_assembly(name, sample_folder, R1, R2, threads):
	## check if previously assembled and succeeded
	filename_stamp = sample_folder + '/.success'
	if os.path.isfile(filename_stamp):
		stamp =	functions.read_time_stamp(filename_stamp)
		print ("\tA previous command generated results on: ", stamp)
	else:
		# Call spades_assembler
		spades_assembler.run_module_SPADES(name, sample_folder, R1, R2, threads)

