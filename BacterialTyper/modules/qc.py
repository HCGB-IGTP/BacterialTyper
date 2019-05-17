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

## import my modules
from BacterialTyper.modules import sample_prepare
from BacterialTyper import fastqc_caller
from BacterialTyper import multiQC_report
from BacterialTyper import BUSCO_caller
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
	start_time_partial = start_time_total
	
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
	start_time_partial = functions.timestamp(start_time_partial)

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
	functions.boxymcboxface("BUSCO Assembly Quality check")
	print ("--------- Starting Process ---------")
	functions.print_time()

	## absolute path for in & out
	input_dir = os.path.abspath(options.input)
	outdir = os.path.abspath(options.output_folder)
	database_folder = os.path.abspath(options.database)

	## get files
	print ("+ Retrieve all genomes assembled...")
	my_assembly_list = functions.retrieve_files(input_dir, '_chromosome.fna')

	## generate output folder
	functions.create_folder(outdir)
	
	## Check each 
	BUSCO_database = database_folder + '/BUSCO'
	dataFrame_results = BUSCO_call(options.BUSCO_dbs, my_assembly_list, BUSCO_database, outdir, options.threads)
	
	print ("+ Quality control of all samples finished: ")
	start_time_partial = functions.timestamp(start_time_partial_BUSCO)
	
	## generate plots
	print ("+ Generate summarizing plots...")
	BUSCO_plots(dataFrame_results, outdir, options.threads)	

################################################
def BUSCO_call(datasets, list_scaffolds, database_folder, output, threads):

	## get datasets
	print ("+ Check folder provided as database for available BUSCO datasets...")
	BUSCO_datasets = BUSCO_caller.BUSCO_retrieve_sets(datasets, database_folder)
	
	## optimize threads
	threads_module = functions.optimize_threads(threads, len(list_scaffolds))

	## BUSCO needs to chdir to output folder
	path_here = os.getcwd()
	
	print ("+ Checking quality for each sample retrieved...")
	
	# We can use a with statement to ensure threads are cleaned up promptly
	with concurrent.futures.ThreadPoolExecutor(max_workers=1) as executor: ## need to do 1 by one as there is a problem with the working directory
		for DataSet in BUSCO_datasets:
			## send for each sample
			commandsSent = { executor.submit( BUSCO_runner, DataSet, sample, BUSCO_datasets[DataSet], output, threads): sample for sample in list_scaffolds }

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
	
	short_summary = pd.DataFrame(columns=('Sample', 'Dataset', 'Summary File', 'Folder'))

	## generate report
	for DataSet in BUSCO_datasets:
		for sample in list_scaffolds:
			sample_name = os.path.basename(sample).split('_chromosome.fna')[0]
			my_BUSCO_results_folder = output + '/' + sample_name + '/run_' + DataSet
			my_short_tsv = 	my_BUSCO_results_folder + '/short_summary_' + DataSet + '.txt'
			
			if os.path.isfile(my_short_tsv):
				short_summary.loc[len(short_summary)] = (sample_name, DataSet, my_short_tsv, my_BUSCO_results_folder)
			
	return (short_summary)

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

################################################
def BUSCO_runner(DataSet, sample, dataset_path, output, threads_module):
	file_name = os.path.basename(sample)
	sample_name = file_name.split('_chromosome.fna')[0]
	sample_name_dir = functions.create_subfolder(sample_name, output)

	## run busco	
	code = BUSCO_caller.BUSCO_run( dataset_path, sample, threads_module, sample_name_dir, DataSet)
	
	if (code == 'FAIL'):
		BUSCO_caller.BUSCO_run( dataset_path, sample, threads_module, sample_name_dir, DataSet)
	else:
		return ()


	


