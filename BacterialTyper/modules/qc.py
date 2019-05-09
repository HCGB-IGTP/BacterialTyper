#!/usr/bin/env python3
'''
This code calls 
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
from BacterialTyper.modules import sample_prepare
from BacterialTyper import fastqc_caller
from BacterialTyper import multiQC_report
from BacterialTyper import functions
from BacterialTyper import config

################################################
def fastqc(options):
	
	## init time
	start_time_total = time.time()
	
	## debugging messages
	global Debug
	if (options.debug):
		Debug = True
	else:
		Debug = False	
	
	## 
	functions.pipeline_header()
	functions.boxymcboxface("Quality check")
	print ("--------- Starting Process ---------")
	functions.print_time()

	## absolute path for in & out
	input_dir = os.path.abspath(options.input)
	outdir = os.path.abspath(options.output_folder)

	## get files
	pd_samples_retrieved = sample_prepare.get_files(options, input_dir)
	## generate output folder
	functions.create_folder(outdir)

	## optimize threads
	## threads_module = functions.optimize_threads(options.threads, pd_samples_retrieved.index.size) ## no threads implementation here. 1 CPU each
	threads_module = 1
	
	print ("+ Checking quality for each sample retrieved...")	
	
	if (options.pair):
		# We can use a with statement to ensure threads are cleaned up promptly
		with concurrent.futures.ThreadPoolExecutor(max_workers=int(options.threads)) as executor:
			commandsSent = { executor.submit(fastqc_caller.run_module_fastqc, outdir, row['R1'], row['R2'], row['samples']): index for index, row in pd_samples_retrieved.iterrows() }
			for cmd2 in concurrent.futures.as_completed(commandsSent):
				details = commandsSent[cmd2]
				try:
					data = cmd2.result()
				except Exception as exc:
					print ('***ERROR:')
					print (cmd2)
					print('%r generated an exception: %s' % (details, exc))
	
	else:
		print ('+ No implementation yet. Sorry.')
		exit()

	print ("+ FASTQC for samples has finished...")

	## functions.timestamp
	start_time_partial = functions.timestamp(start_time_total)

	if (options.skip_report):
		print ("+ No report generation...")
	else:
		## get subdirs generated and call multiQC report module
		givenList = []
		print ("+ Detail information for each sample could be identified in separate folders:")
		for index, row in pd_samples_retrieved.iterrows():
			fold_sample = outdir + '/' + row['samples']
			givenList.append(fold_sample)
			print ('\t- %s' %fold_sample)	
		outdir_report = functions.create_subfolder("report", outdir)
		multiQC_report.multiQC_module_call(givenList, "FASTQC", outdir_report)

		print ('\n+ A summary HTML report of each sample is generated in folder: %s' %outdir_report)

	print ("\n*************** Finish *******************")
	start_time_partial = functions.timestamp(start_time_total)

	print ("+ Exiting qc module.")
	exit()
	
################################################
def assembly_check(options):
	print ("Hi there! lets call BUSCO\n")
	print ("To be implemented...\nNot available yet...Sorry!\n")

	# spades_assembler.stats(fasta_file)
