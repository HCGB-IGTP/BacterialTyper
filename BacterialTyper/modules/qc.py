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
import pandas as pd
import shutil
from termcolor import colored

## import my modules
from BacterialTyper.modules import sample_prepare
from BacterialTyper import fastqc_caller
from BacterialTyper import sampleParser
from BacterialTyper import multiQC_report
from BacterialTyper import BUSCO_caller
from BacterialTyper import functions
from BacterialTyper import config
from BacterialTyper.modules import info

################################################
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
	
	## set main header
	functions.pipeline_header()
	functions.boxymcboxface("Quality check")
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
	
	### option
	if (options.raw_reads):
		fastqc(input_dir, outdir, options, start_time_total)
	elif (options.assembly):
		BUSCO_check(input_dir, outdir, options, start_time_total, "genome")
	elif (options.annotation):
		BUSCO_check(input_dir, outdir, options, start_time_total, "proteins")

################################################
def fastqc(input_dir, outdir, options, start_time_total):
	
	functions.boxymcboxface("FASTQC Quality check for samples")
	
	## get files
	pd_samples_retrieved = sample_prepare.get_files(options, input_dir, "fastq", "raw")
	
	## debug message
	if (Debug):
		print (colored("\n**DEBUG: pd_samples_retrieve **", 'yellow'))
		print (pd_samples_retrieved)
		print ("\n")

	## generate output folder, if necessary
	print ("\n+ Create output folder(s):")
	outdir_dict = functions.outdir_project(outdir, options.project, pd_samples_retrieved, "fastqc")
	
	if not options.project:
		functions.create_folder(outdir)

	## optimize threads
	## threads_module = functions.optimize_threads(options.threads, pd_samples_retrieved.index.size) ## no threads implementation here. 1 CPU each
	threads_module = 1
	
	print ("+ Checking quality for each sample retrieved...")
	start_time_partial = start_time_total
	
	# Group dataframe by sample name
	sample_frame = pd_samples_retrieved.groupby(["name"])

	# We can use a with statement to ensure threads are cleaned up promptly
	with concurrent.futures.ThreadPoolExecutor(max_workers=int(options.threads)) as executor:
		commandsSent = { executor.submit(fastqc_caller.run_module_fastqc, outdir_dict[name], sorted( cluster["sample"].tolist() ), name): name for name, cluster in sample_frame }
		
		for cmd2 in concurrent.futures.as_completed(commandsSent):
			details = commandsSent[cmd2]
			try:
				data = cmd2.result()
			except Exception as exc:
				print ('***ERROR:')
				print (cmd2)
				print('%r generated an exception: %s' % (details, exc))

	print ("+ FASTQC for samples has finished...")	
	
	## functions.timestamp
	start_time_partial = functions.timestamp(start_time_partial)

	if (options.skip_report):
		print ("+ No report generation...")
	else:
		print ("\n+ Generating a report using MultiQC module.")
		outdir_report = functions.create_subfolder("report", outdir)

		## get subdirs generated and call multiQC report module
		givenList = []
		print ("+ Detail information for each sample could be identified in separate folders:")
		
		## call multiQC report module
		givenList = [ v for v in outdir_dict.values() ]
		my_outdir_list = set(givenList)

		## debug message
		if (Debug):
			print (colored("\n**DEBUG: my_outdir_list for multiqc report **", 'yellow'))
			print (my_outdir_list)
			print ("\n")
		
		fastqc_report = functions.create_subfolder("FASTQC", outdir_report)
		multiQC_report.multiQC_module_call(my_outdir_list, "FASTQC", fastqc_report,"")
		print ('\n+ A summary HTML report of each sample is generated in folder: %s' %fastqc_report)

	print ("\n*************** Finish *******************")
	start_time_partial = functions.timestamp(start_time_total)

	print ("+ Exiting qc module.")
	exit()

################################################
def BUSCO_check(input_dir, outdir, options, start_time_total, mode):

	functions.boxymcboxface("BUSCO Assembly Quality check")

	## absolute path for in & out
	database_folder = os.path.abspath(options.database)

	## get files
	pd_samples_retrieved = sample_prepare.get_files(options, input_dir, "assemble", "fasta")

	## debug message
	if (Debug):
		print (colored("**DEBUG: df_samples_busco **", 'yellow'))
		print (df_samples_busco)

	exit()

	## get dir for each sample
	BUSCO_outdir_dict = functions.outdir_project(outdir, options.project, pd_samples_retrieved, "busco_qc")
	
	## initiate dataframe
	name_columns = ("name", "assembly", "busco_qc")
	df_samples_busco = pd.DataFrame(columns=name_columns)
	for f in outdir_dict:
		my_assembly_file = spades_assembler.get_files(outdir_dict[f])
		df_samples_busco.loc[len(df_samples_busco)] = [(f, my_assembly_file, BUSCO_outdir_dict[f])]
		
	## debug message
	if (Debug):
		print (colored("**DEBUG: df_samples_busco **", 'yellow'))
		print (df_samples_busco)

	exit()

	## Check each using BUSCO
	functions.boxymcboxface("BUSCO Assembly Quality check")
	database_folder = os.path.abspath(options.database)
	dataFrame_results = BUSCO_call(options.BUSCO_dbs, df_samples_busco, options.threads, "genome")
	
	## functions.timestamp
	print ("+ Quality control of all samples finished: ")
	start_time_partial = functions.timestamp(start_time_partial_BUSCO)
	
	## multiqc report plot
	if (options.skip_report):
		print ("+ No report generation...")
	else:
		## generate plots
		print ("+ Generate summarizing plots...")
		BUSCO_plots(dataFrame_results, outdir, options.threads)	

################################################
def BUSCO_call(datasets, pd_samples, database_folder, output, threads, mode):

	## mode= proteins|genome

	## get datasets
	print ("+ Check folder provided as database for available BUSCO datasets...")
	BUSCO_datasets = BUSCO_caller.BUSCO_retrieve_sets(datasets, database_folder)
	
	## BUSCO needs to chdir to output folder
	path_here = os.getcwd()
	
	print ("+ Checking quality for each sample retrieved...")
	
	## optimize threads: No need to optimize. There is a problem with the working dir of BUSCO and we need to change every time
	## threads_module = functions.optimize_threads(threads, len(list_scaffolds))

	# We can use a with statement to ensure threads are cleaned up promptly
	with concurrent.futures.ThreadPoolExecutor(max_workers=1) as executor: ## need to do 1 by one as there is a problem with the working directory
		for DataSet in BUSCO_datasets:
			## send for each sample
			commandsSent = { executor.submit( BUSCO_runner, name, DataSet, row['sample'], BUSCO_datasets[DataSet], output, threads, mode): name for name, row in pd_samples }

			for cmd2 in concurrent.futures.as_completed(commandsSent):
				details = commandsSent[cmd2]
				try:
					data = cmd2.result()
				except Exception as exc:
					print ('***ERROR:')
					print (cmd2)
					print('%r generated an exception: %s' % (details, exc))
			
			print ("+ Jobs finished for dataset %s\n+ Collecting information..." %DataSet)

	print ("Finish here...")
	os.chdir(path_here)
	
	## init dataframe
	short_summary = pd.DataFrame(columns=('Sample', 'Dataset', 'Summary File', 'Folder'))

	## generate report
	for DataSet in BUSCO_datasets:
		for sample in list_scaffolds:
			file_name = os.path.basename(sample)
			sample_name =""
			if file_name.endswith('_chromosome.fna'): ## mode genome
				sample_name = file_name.split('_chromosome.fna')[0]
			else: 
				sample_name = file_name.split('.faa')[0] ## mode proteins
		
			my_BUSCO_results_folder = output + '/' + sample_name + '/run_' + DataSet
			my_short_tsv = 	my_BUSCO_results_folder + '/short_summary_' + DataSet + '.txt'
			
			if os.path.isfile(my_short_tsv):
				short_summary.loc[len(short_summary)] = (sample_name, DataSet, my_short_tsv, my_BUSCO_results_folder)
			
	return (short_summary)

################################################
def BUSCO_runner(sample_name, DataSet, sample, dataset_path, output, threads_module, mode):
	#file_name = os.path.basename(sample)
	sample_name_dir = functions.create_subfolder(sample_name, output)

	## run busco	
	code = BUSCO_caller.BUSCO_run( dataset_path, sample, threads_module, sample_name_dir, DataSet, mode)
	
	if (code == 'FAIL'):
		BUSCO_caller.BUSCO_run( dataset_path, sample, threads_module, sample_name_dir, DataSet, mode)
	else:
		return ()

################################################
def BUSCO_plots(dataFrame_results, outdir, threads):
	
	list_datasets = set(dataFrame_results['Dataset'].tolist())
	list_samples = set(dataFrame_results['Sample'].tolist())

	plot_folder = functions.create_subfolder('BUSCO_plots', outdir)
	outdir_busco_plot = []
	
	print ("+ Get results for all samples summarized by dataset:")
	for dataset in list_datasets:
		print ("\t+ Get results for: ", dataset)
		plot_folder_dataset = functions.create_subfolder(dataset, plot_folder)
		outdir_busco_plot.append(plot_folder_dataset)
	
		for index, row in dataFrame_results.iterrows():
			if (dataset == row['Dataset']):
				shutil.copy(row['Summary File'], plot_folder_dataset + '/short_summary_' + dataset + '_' + row['Sample'] + '.txt')
		
	print ("+ Get results for summarized by sample:")
	for sample in list_samples:
		print ("\t+ Get results for: ", sample)
		plot_folder_sample = functions.create_subfolder(sample, plot_folder)
		outdir_busco_plot.append(plot_folder_sample)

		for index, row in dataFrame_results.iterrows():
			if (sample == row['Sample']):
				shutil.copy(row['Summary File'], plot_folder_sample + '/short_summary_' + row['Dataset'] + '_' + sample + '.txt')
	
	print ("+ Generate plots for each subset")
	
	# We can use a with statement to ensure threads are cleaned up promptly
	with concurrent.futures.ThreadPoolExecutor(max_workers=threads) as executor: ## need to do 1 by one as there is a problem with the working directory
		## send for each sample
		commandsSent = { executor.submit( BUSCO_caller.BUSCO_plot , plot): plot for plot in outdir_busco_plot }

		for cmd2 in concurrent.futures.as_completed(commandsSent):
			details = commandsSent[cmd2]
			try:
				data = cmd2.result()
			except Exception as exc:
				print ('***ERROR:')
				print (cmd2)
				print('%r generated an exception: %s' % (details, exc))
			
	print ("+ All plots generated...")
	print ("+ Check results under folders in : ", plot_folder)		

