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
import concurrent.futures
from termcolor import colored

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
		
	##################################
	### show help messages if desired	
	##################################
	if (options.help_format):
		## help_format option
		sampleParser.help_format()
		exit()
	elif (options.help_BUSCO):
		## information for BUSCO
		BUSCO_caller.print_help_BUSCO()
		exit()
	elif (options.help_project):
		## information for project
		help.project_help()
		exit()
	elif (options.help_multiqc):
		## information for Multiqc
		multiQC_report.multiqc_help()
		exit()
	
	### set as default paired_end mode
	if (options.single_end):
		options.pair = False
	else:
		options.pair = True

	## message header
	functions.pipeline_header()
	functions.boxymcboxface("Assembly module")
	print ("--------- Starting Process ---------")
	functions.print_time()
	
	## absolute path for in & out
	input_dir = os.path.abspath(options.input)
	outdir=""

	## set mode: project/detached
	if (options.project):
		outdir = input_dir		
	elif (options.detached):
		outdir = os.path.abspath(options.output_folder)

	## get files
	pd_samples_retrieved = sample_prepare.get_files(options, input_dir, "fastq", "_trim_")
	
	## debug message
	if (Debug):
		print (colored("**DEBUG: pd_samples_retrieve **", 'yellow'))
		print (pd_samples_retrieved)
	
	## generate output folder, if necessary
	print ("\n+ Create output folder(s):")
	outdir_dict = functions.outdir_project(outdir, options.project, pd_samples_retrieved, "assemble")
	
	if not options.project:
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

	# Group dataframe by sample name
	sample_frame = pd_samples_retrieved.groupby(["name"])

	# We can use a with statement to ensure threads are cleaned up promptly
	with concurrent.futures.ThreadPoolExecutor(max_workers=int(workers)) as executor:
		## send for each sample
		commandsSent = { executor.submit( check_sample_assembly, name, outdir_dict[name],  sorted(cluster["sample"].tolist()), threads_module): name for name, cluster in sample_frame }

		for cmd2 in concurrent.futures.as_completed(commandsSent):
			details = commandsSent[cmd2]
			try:
				data = cmd2.result()
			except Exception as exc:
				print ('***ERROR:')
				print (cmd2)
				print('%r generated an exception: %s' % (details, exc))
		
	## functions.timestamp
	print ("+ Assembly of all samples finished: ")
	start_time_partial = functions.timestamp(start_time_partial_assembly)
	start_time_partial_BUSCO = start_time_partial

	### symbolic links
	print ("+ Retrieve all genomes assembled...")
	
	### BUSCO check assembly
	qc.BUSCO_check(input_dir, outdir, options, start_time_total_BUSCO, "genome")
	
	print ("\n*************** Finish *******************")
	start_time_partial = functions.timestamp(start_time_total)

	print ("+ Exiting Assembly module.")
	exit()

#############################################
def check_sample_assembly(name, sample_folder, files, threads):
	## check if previously assembled and succeeded
	filename_stamp = sample_folder + '/.success'
	if os.path.isfile(filename_stamp):
		stamp =	functions.read_time_stamp(filename_stamp)
		print (colored("\tA previous command generated results on: %s" %stamp, 'yellow'))
	else:
	
		## debug message
		if (Debug):
			print (colored("**DEBUG: spades_assembler.run_module_SPADES call**", 'yellow'))
			print ("spades_assembler.run_module_SPADES(name, sample_folder, files[0], files[1], threads)")
			print ("spades_assembler.run_module_SPADES " + name + "\t" + sample_folder + "\t" + files[0] + "\t" + files[1] + "\t" + str(threads) + "\n")
	
		# Call spades_assembler
		spades_assembler.run_module_SPADES(name, sample_folder, files[0], files[1], threads)

