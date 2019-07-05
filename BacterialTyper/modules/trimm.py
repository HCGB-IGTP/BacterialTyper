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
import time
from io import open
import shutil
import concurrent.futures
from termcolor import colored

## import my modules
from BacterialTyper import trimmomatic_call
from BacterialTyper import multiQC_report
from BacterialTyper.modules import sample_prepare
from BacterialTyper import sampleParser
from BacterialTyper import functions
from BacterialTyper import config

##############################################
def run(options):

	## init time
	start_time_total = time.time()

	##################################
	### show help messages if desired	
	##################################
	if (options.help_format):
		## help_format option
		sampleParser.help_format()
		exit()
	elif (options.help_trimm_adapters):
		## help on trimm adapters
		trimmomatic_call.print_help_adapters()
		exit()
	elif (options.help_project):
		## information for project
		help.project_help()
		exit()
	elif (options.help_multiqc):
		## information for Multiqc
		multiQC_report.multiqc_help()
		exit()
		
	## debugging messages
	global Debug
	if (options.debug):
		Debug = True
	else:
		Debug = False
		
	### set as default paired_end mode
	if (options.single_end):
		options.pair = False
	else:
		options.pair = True
	
	functions.pipeline_header()
	functions.boxymcboxface("Trimming samples")
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
	pd_samples_retrieved = sample_prepare.get_files(options, input_dir, "fastq", ("fastq", "fq", "fastq.gz", "fq.gz"))
	
	## debug message
	if (Debug):
		print (colored("**DEBUG: pd_samples_retrieve **", 'yellow'))
		print (pd_samples_retrieved)

	## generate output folder, if necessary
	print ("\n+ Create output folder(s):")
	if not options.project:
		functions.create_folder(outdir)
	## for samples
	outdir_dict = functions.outdir_project(outdir, options.project, pd_samples_retrieved, "trimm")
	
	## optimize threads
	## workers:
	name_list = set(pd_samples_retrieved["name"].tolist())
	max_workers_int = len(name_list)
	
	## number_samples = pd_samples_retrieved.index.size => Number samples
	threads_module = functions.optimize_threads(options.threads, max_workers_int) ## fix threads_module optimization

	## debug message
	if (Debug):
		print (colored("**DEBUG: options.threads " +  str(options.threads) + " **", 'yellow'))
		print (colored("**DEBUG: max workers " +  str(max_workers_int) + " **", 'yellow'))
		print (colored("**DEBUG: cpu/process " +  str(threads_module ) + " **", 'yellow'))

	print ("+ Trimming adapters for each sample retrieved...")	
	
	# Group dataframe by sample name
	sample_frame = pd_samples_retrieved.groupby(["name"])
	
	## send for each sample
	with concurrent.futures.ThreadPoolExecutor(max_workers=int(max_workers_int)) as executor:
		commandsSent = { executor.submit(trimmo_caller, sorted(cluster["sample"].tolist()), outdir_dict[name], name, threads_module, Debug, options.adapters): name for name, cluster in sample_frame }

		for cmd2 in concurrent.futures.as_completed(commandsSent):
			details = commandsSent[cmd2]
			try:
				data = cmd2.result()
			except Exception as exc:
				print ('***ERROR:')
				print (cmd2)
				print('%r generated an exception: %s' % (details, exc))

	print ("\n\n+ Trimming samples has finished...")
	## functions.timestamp
	start_time_partial = functions.timestamp(start_time_total)

	## get files generated and generate symbolic link
	if not options.project:
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
		outdir_report = functions.create_subfolder("report", outdir)
	
		## call multiQC report module
		givenList = [ v for v in outdir_dict.values() ]
		my_outdir_list = set(givenList)
		
		## debug message
		if (Debug):
			print (colored("\n**DEBUG: my_outdir_list for multiqc report **", 'yellow'))
			print (my_outdir_list)
			print ("\n")

		trimm_report = functions.create_subfolder("trimm", outdir_report)
		multiQC_report.multiQC_module_call(my_outdir_list, "Trimmomatic", trimm_report,"")
		print ('\n+ A summary HTML report of each sample is generated in folder: %s' %trimm_report)
		
	print ("\n*************** Finish *******************")
	start_time_partial = functions.timestamp(start_time_total)
	print ("\n+ Exiting trimm module.")
	exit()
	

#############################################
def trimmo_caller(list_reads, sample_folder, name, threads, Debug, adapters):
	## check if previously assembled and succeeded
	filename_stamp = sample_folder + '/.success'
	if os.path.isfile(filename_stamp):
		stamp =	functions.read_time_stamp(filename_stamp)
		print (colored("\tA previous command generated results on: %s" %stamp, 'yellow'))
	else:
		# Call trimmomatic
		trimmomatic_call.trimmo_module(list_reads, sample_folder, name, threads, Debug, adapters)



	
