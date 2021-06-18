#!/usr/bin/env python3
##########################################################
## Jose F. Sanchez										##
## Copyright (C) 2019 Lauro Sumoy Lab, IGTP, Spain		##
##########################################################
"""
Trimms sequence adapters within fastq reads.
"""
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
from BacterialTyper.scripts import trimmomatic_call
from BacterialTyper.scripts import multiQC_report
from BacterialTyper.config import set_config
from BacterialTyper.modules import help_info
from BacterialTyper.modules import qc
from BacterialTyper.data import data_files
from BacterialTyper import __version__ as pipeline_version

import HCGB
from HCGB import sampleParser
import HCGB.functions.aesthetics_functions as HCGB_aes
import HCGB.functions.time_functions as HCGB_time
import HCGB.functions.main_functions as HCGB_main
import HCGB.functions.files_functions as HCGB_files
import HCGB.functions.info_functions as HCGB_info

##############################################
def run(options):

	## init time
	start_time_total = time.time()

	##################################
	### show help messages if desired	
	##################################
	if (options.help_format):
		## help_format option
		help_info.help_fastq_format()
		exit()
	elif (options.help_trimm_adapters):
		## help on trimm adapters
		trimmomatic_call.print_help_adapters()
		exit()
	elif (options.help_project):
		## information for project
		help_info.project_help()
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
	
	HCGB_aes.pipeline_header("BacterialTyper", ver=pipeline_version)
	HCGB_aes.boxymcboxface("Trimming samples")
	print ("--------- Starting Process ---------")
	HCGB_time.print_time()

	## absolute path for in & out
	input_dir = os.path.abspath(options.input)
	options.input = input_dir
	outdir=""

	## Default parameters:
	trimmomatic_params = {
		"ILLUMINACLIP": options.ILLUMINACLIP,
		"LEADING": str(options.LEADING),
		"TRAILING": str(options.TRAILING),
		"SLIDINGWINDOW": options.SLIDINGWINDOW,
		"MINLEN": str(options.MINLEN)
		}

	## Project mode as default
	if (options.detached):
		options.project = False
		outdir = os.path.abspath(options.output_folder)
	else:
		options.project = True
		outdir = input_dir	
	
	##
	options.output_folder = outdir
	
	## get files
	pd_samples_retrieved = sampleParser.files.get_files(options, input_dir, "fastq", ("fastq", "fq", "fastq.gz", "fq.gz"), options.debug)
	
	## debug message
	if (Debug):
		HCGB_aes.debug_message("pd_samples_retrieved", 'yellow')
		HCGB_main.print_all_pandaDF(pd_samples_retrieved)

	## generate output folder, if necessary
	print ("\n+ Create output folder(s):")
	if not options.project:
		HCGB_files.create_folder(outdir)
	## for samples
	outdir_dict = HCGB_files.outdir_project(outdir, options.project, pd_samples_retrieved, "trimm", options.debug)
	
	## optimize threads
	name_list = set(pd_samples_retrieved["name"].tolist())
	threads_job = HCGB_main.optimize_threads(options.threads, len(name_list)) ## threads optimization
	max_workers_int = int(options.threads/threads_job)

	## debug message
	if (Debug):
		print (colored("**DEBUG: options.threads " +  str(options.threads) + " **", 'yellow'))
		print (colored("**DEBUG: max_workers " +  str(max_workers_int) + " **", 'yellow'))
		print (colored("**DEBUG: cpu_here " +  str(threads_job) + " **", 'yellow'))

	print ("+ Trimming adapters for each sample retrieved...")	
	
	# Group dataframe by sample name
	sample_frame = pd_samples_retrieved.groupby(["name"])
	
	# Trimming adapters
	if (options.adapters):
		# Adapter file provided
		options.adapters = os.path.abspath(options.adapters)
		print("\t- Adapters file provided...")
	else:
		# Get default adpaters file
		print("\t- Default Trimmomatic adapters (v0.39) will be used...")
		options.adapters = data_files.data_list("available_Trimmomatic_adapters")
		
	## send for each sample
	with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers_int) as executor:
		commandsSent = { executor.submit(trimmo_caller, sorted(cluster["sample"].tolist()), outdir_dict[name], name, threads_job, Debug, trimmomatic_params, options.adapters): name for name, cluster in sample_frame }

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
	start_time_partial = HCGB_time.timestamp(start_time_total)

	## get files generated and generate symbolic link
	if not options.project:
		dir_symlinks = HCGB_files.create_subfolder('link_files', outdir)
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
	
		HCGB_files.get_symbolic_link(files2symbolic, dir_symlinks)

	if (options.skip_report):
		print ("+ No report generation...")
	else:
		print ("\n+ Generating a report using MultiQC module.")
		outdir_report = HCGB_files.create_subfolder("report", outdir)
	
		## call multiQC report module
		givenList = [ v for v in outdir_dict.values() ]
		my_outdir_list = set(givenList)
		
		## debug message
		if (Debug):
			HCGB_aes.debug_message("my_outdir_list for multiqc report", "yellow")
			print (my_outdir_list)
			print ("\n")

		trimm_report = HCGB_files.create_subfolder("trimm", outdir_report)
		multiQC_report.multiQC_module_call(my_outdir_list, "Trimmomatic", trimm_report,"")
		print ('\n+ A summary HTML report of each sample is generated in folder: %s' %trimm_report)
		
		## create fastqc for trimmed reads
		pd_samples_retrieved_trimmed = sampleParser.files.get_files(options, input_dir, "trim", ['_trim'], options.debug)
		qc.fastqc(pd_samples_retrieved_trimmed, outdir, options, start_time_partial, "trimmed", Debug)
		
	print ("\n*************** Finish *******************")
	start_time_partial = HCGB_time.timestamp(start_time_total)
	
	################################################
	## dump information and parameters
	################################################
	## samples information dictionary
	samples_info = {}
	samples_frame = pd_samples_retrieved_trimmed.groupby('new_name')
	for name, grouped in samples_frame:
		samples_info[name] = grouped['sample'].to_list()
	
	## trimmomatic params
	del options.MINLEN
	del options.ILLUMINACLIP
	del options.LEADING
	del options.TRAILING
	del options.SLIDINGWINDOW
	
	info_dir = HCGB_files.create_subfolder("info", outdir)
	print("+ Dumping information and parameters")
	runInfo = { "module":"trimm", "time":time.time(),
				"BacterialTyper version":pipeline_version,
				'sample_info': samples_info,
				'trimmomatic_params': trimmomatic_params }
	
	HCGB_info.dump_info_run(info_dir, 'trimm', options, runInfo, options.debug)

	print ("\n+ Exiting trimm module.")
	return()
	
#############################################
def trimmo_caller(list_reads, sample_folder, name, threads, Debug, trimmomatic_params, adapters):
	## check if previously assembled and succeeded
	filename_stamp = sample_folder + '/.success'
	if os.path.isfile(filename_stamp):
		stamp =	HCGB_time.read_time_stamp(filename_stamp)
		print (colored("\tA previous command generated results on: %s [%s]" %(stamp, name), 'yellow'))
	else:
		# Call trimmomatic
		trimmomatic_call.trimmo_module(list_reads, sample_folder, name, threads, Debug, adapters, trimmomatic_params)
