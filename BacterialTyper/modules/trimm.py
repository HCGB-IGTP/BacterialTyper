#!/usr/bin/env python3
'''
This code calls Trimmomatic for the trimming of sequence adapter within fastq reads.
Jose F. Sanchez
Copyright (C) 2019 Lauro Sumoy Lab, IGTP, Spain
'''
## import useful modules
import os
import sys
import re
from io import open
import shutil
import concurrent.futures

## import my modules
from BacterialTyper import trimmomatic_call
from BacterialTyper import multiQC_report
from BacterialTyper.modules import sample_prepare
from BacterialTyper import functions
from BacterialTyper import config

##############################################
def run(options):

	## debugging messages
	global Debug
	if (options.debug):
		Debug = True
	else:
		Debug = False
	
	functions.pipeline_header()
	functions.boxymcboxface("Trimming samples")
	
	## absolute path for in & out
	input_dir = os.path.abspath(options.input)
	outdir = os.path.abspath(options.output_folder)

	## get files
	pd_samples_retrieved = sample_prepare.get_files(options, input_dir)
	
	## debug message
	if (Debug):
		print (colored("**DEBUG: pd_samples_retrieve **", 'yellow'))
		print (pd_samples_retrieved)

	## generate output folder
	functions.create_folder(outdir)

	## optimize threads
	## number_samples = pd_samples_retrieved.index.size => Number samples
	threads_module = functions.optimize_threads(options.threads, pd_samples_retrieved.index.size)
	
	print ("+ Trimming adapters for each sample retrieved...")	
	
	if (options.pair):
		# We can use a with statement to ensure threads are cleaned up promptly
		with concurrent.futures.ThreadPoolExecutor(max_workers=int(options.threads)) as executor:
			commandsSent = { executor.submit(trimmomatic_call.trimmo_module, row['R1'], row['R2'], outdir, row['samples'], threads_module, Debug): index for index, row in pd_samples_retrieved.iterrows() }
			for cmd2 in concurrent.futures.as_completed(commandsSent):
				details = commandsSent[cmd2]
				try:
					data = cmd2.result()
				except Exception as exc:
					print ('***ERROR:')
					print (cmd2)
					print('%r generated an exception: %s' % (details, exc))
	
	else:
		# We can use a with statement to ensure threads are cleaned up promptly
		with concurrent.futures.ThreadPoolExecutor(max_workers=num_threads) as executor:
			commandsSent = { executor.submit(trimmomatic_call.trimmo_module, row['read'], 'na', outdir, row['samples'], threads_module, Debug): index for index, row in pd_samples_retrieved.iterrows() }
			for cmd2 in concurrent.futures.as_completed(commandsSent):
				details = commandsSent[cmd2]
				try:
					data = cmd2.result()
				except Exception as exc:
					print ('***ERROR:')
					print (cmd2)
					print('%r generated an exception: %s' % (details, exc))


	print ("\n\n+ Trimming samples has finished...")
	
	## get files generated and generate symbolic link
	dir_symlinks = functions.create_subfolder('link_files', outdir)
	files2symbolic = []
	folders = os.listdir(outdir)
	## debug message
	if (Debug):
		print (colored("**DEBUG: generate symbolic links for each file in " + dir_symlinks + "**", 'yellow'))
		
	for fold in folders:
		if fold.endswith(".log"):
			continue
		else:
			this_folder = outdir + '/' + fold
			subfiles = os.listdir(this_folder)
			for files in subfiles:
				files_search = re.search(r".*trim_R\d{1}.*", files) ## only paired-end. Todo: single end
				if files_search:
					files2symbolic.append(this_folder + '/' + files)
	
	functions.get_symbolic_link(files2symbolic, dir_symlinks)
					
	if (options.skip_report):
		print ("+ No report generation...")
	else:
		print ("\n+ Generating a report using MultiQC module.")
	
		## call multiQC report module
		givenList = []
		givenList.append(outdir)
		my_outdir_list = set(givenList)
		outdir_report = functions.create_subfolder("report", outdir)
		multiQC_report.multiQC_module_call(my_outdir_list, "Trimmomatic", outdir_report)

	print ("\n+ Exiting trimm module.")
	exit()
