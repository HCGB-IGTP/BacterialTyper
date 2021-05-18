#!/usr/bin/env python3
##########################################################
## Jose F. Sanchez										##
## Copyright (C) 2019 Lauro Sumoy Lab, IGTP, Spain		##
##########################################################
'''
Generates sample quality control at raw reads, assembly or annotation level.
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

## import my modules, scripts, config
from BacterialTyper.scripts import fastqc_caller
from BacterialTyper.scripts import multiQC_report
from BacterialTyper.scripts import BUSCO_caller
from BacterialTyper.config import set_config
from BacterialTyper.modules import help_info

import HCGB
from HCGB import sampleParser
import HCGB.functions.aesthetics_functions as HCGB_aes
import HCGB.functions.time_functions as HCGB_time
import HCGB.functions.main_functions as HCGB_main
import HCGB.functions.files_functions as HCGB_files

################################################
def run_QC(options):

	## init time
	start_time_total = time.time()

	##################################
	### show help messages if desired	
	##################################
	if (options.help_format):
		## help_format option
		help_info.help_fastq_format()
		exit()
	elif (options.help_BUSCO):
		## information for BUSCO
		BUSCO_caller.print_help_BUSCO()
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
	
	## set main header
	HCGB_aes.pipeline_header("BacterialTyper")
	HCGB_aes.boxymcboxface("Quality check")
	print ("--------- Starting Process ---------")
	HCGB_time.print_time()

	## absolute path for in & out
	input_dir = os.path.abspath(options.input)
	outdir=""

	## Project mode as default
	if (options.detached):
		options.project = False
		outdir = os.path.abspath(options.output_folder)
	else:
		options.project = True
		outdir = input_dir		
	
	### option
	if (options.raw_reads):
		## get files
		pd_samples_retrieved = sampleParser.files.get_files(options, input_dir, "fastq", ("fastq", "fq", "fastq.gz", "fq.gz"), options.debug)
		fastqc(pd_samples_retrieved, outdir, options, start_time_total, "raw", Debug)
	elif (options.trim_reads):
		## get files
		pd_samples_retrieved = sampleParser.files.get_files(options, input_dir, "trim", ['_trim'], options.debug)
		fastqc(pd_samples_retrieved, outdir, options, start_time_total, "trimmed", Debug)
	elif (options.assembly):
		BUSCO_check(input_dir, outdir, options, start_time_total, "genome")
	elif (options.annotation):
		BUSCO_check(input_dir, outdir, options, start_time_total, "proteins")
		
	return()

################################################
def fastqc(pd_samples_retrieved, outdir, options, start_time_total, name_analysis, Debug):
	
	HCGB_aes.boxymcboxface("FASTQC Quality check for samples")
	
	## debug message
	if (Debug):
		print (colored("\n**DEBUG: pd_samples_retrieve **", 'yellow'))
		print (pd_samples_retrieved)
		print ("\n")

	## generate output folder, if necessary
	print ("\n+ Create output folder(s):")
	
	## if not project, outdir contains the dir to put output
	## in this case, in some other cases might not occur	
	if not options.project:
		functions.create_folder(outdir)
	outdir_dict = HCGB_files.outdir_project(outdir, options.project, pd_samples_retrieved, "fastqc_" + name_analysis, options.debug)
	
	print ("+ Checking quality for each sample retrieved...")
	start_time_partial = start_time_total
	
	# Group dataframe by sample name
	sample_frame = pd_samples_retrieved.groupby(["name"])

	## optimize threads
	name_list = set(pd_samples_retrieved["name"].tolist())
	threads_job = HCGB_main.optimize_threads(options.threads, len(name_list)) ## threads optimization
	max_workers_int = int(options.threads/threads_job)

	## debug message
	if (Debug):
		HCGB_aes.debug_message("options.threads: " + str(options.threads), "yellow")
		HCGB_aes.debug_message("max_workers: " + str(max_workers_int), "yellow")
		HCGB_aes.debug_message("threads_job: " + str(threads_job), "yellow")

	## send for each sample
	print ("+ Calling fastqc for samples...")	
	with concurrent.futures.ThreadPoolExecutor(max_workers=int(max_workers_int)) as executor:
		commandsSent = { executor.submit(fastqc_caller.run_module_fastqc, outdir_dict[name], sorted( cluster["sample"].tolist() ), name, threads_job): name for name, cluster in sample_frame }
		
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
	start_time_partial = HCGB_time.timestamp(start_time_partial)

	if (options.skip_report):
		print ("+ No report generation...")
	else:
		print ("\n+ Generating a report using MultiQC module.")
		outdir_report = HCGB_files.create_subfolder("report", outdir)

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
		
		fastqc_report = HCGB_files.create_subfolder("FASTQC", outdir_report)
		fastqc_final_report = HCGB_files.create_subfolder(name_analysis, fastqc_report)
		multiQC_report.multiQC_module_call(my_outdir_list, "FASTQC", fastqc_final_report,"")
		print ('\n+ A summary HTML report of each sample is generated in folder: %s' %fastqc_final_report)

	print ("\n*************** Finish *******************")
	start_time_partial = HCGB_time.timestamp(start_time_total)

	print ("+ Exiting qc module.")
	exit()

################################################
def BUSCO_check(input_dir, outdir, options, start_time_total, mode):

	HCGB_aes.boxymcboxface("BUSCO Assembly Quality check")

	## absolute path for in & out
	database_folder = os.path.abspath(options.database)

	## get files and get dir for each sample according to mode
	if mode == 'genome':
		pd_samples_retrieved = sampleParser.files.get_files(options, input_dir, "assembly", ["fna"], options.debug)

		if not options.project:
			outdir = HCGB_files.create_subfolder("assembly_qc", outdir)

		if options.debug:
			print ("** DEBUG: pd_samples_retrieved")
			print (pd_samples_retrieved)
		
		BUSCO_outdir_dict = HCGB_files.outdir_project(outdir, options.project, pd_samples_retrieved, "assemble_qc", options.debug)

	elif mode == 'proteins':
		pd_samples_retrieved = sampleParser.files.get_files(options, outdir, "annot", ["faa"], options.debug) ##

		if not options.project:
			outdir = HCGB_files.create_subfolder("annot_qc", outdir)

		if options.debug:
			print ("** DEBUG: pd_samples_retrieved")
			print (pd_samples_retrieved)
			
		BUSCO_outdir_dict = HCGB_files.outdir_project(outdir, options.project, pd_samples_retrieved, "annot_qc", options.debug)

	## add column to dataframe
	pd_samples_retrieved['busco_folder'] = ""
	for index, row in pd_samples_retrieved.iterrows():
		pd_samples_retrieved.at[index, 'busco_folder'] = BUSCO_outdir_dict[ row['name'] ] 

	## debug message
	if (options.debug):
		HCGB_aes.debug_message("df_samples_busco", 'yellow')
		print (pd_samples_retrieved)
		
		HCGB_aes.debug_message("BUSCO_outdir_dict", 'yellow')
		print (BUSCO_outdir_dict)

	## Check each using BUSCO
	database_folder = os.path.abspath(options.database)
	BUSCO_Database = database_folder + '/BUSCO'
	if not os.path.exists(BUSCO_Database):
		HCGB_files.create_folder(BUSCO_Database)

	## call
	(dataFrame_results, stats_results) = BUSCO_call(options.BUSCO_dbs, pd_samples_retrieved, BUSCO_Database, options.threads, mode)
	
	## debug message
	if (options.debug):
		HCGB_aes.debug_message("dataFrame_results", 'yellow')
		HCGB_main.print_all_pandaDF(dataFrame_results)
	
	## functions.timestamp
	print ("+ Quality control of all samples finished: ")
	start_time_partial = HCGB_time.timestamp(start_time_total)
	
	## multiqc report plot
	if (options.skip_report):
		print ("+ No report generation...")
	else:
		print ("\n+ Generating a report BUSCO plot.")
		outdir_report = HCGB_files.create_subfolder("report", outdir)

		## get subdirs generated and call multiQC report module
		givenList = []
		print ("+ Detail information for each sample could be identified in separate folders.")
		
		## name folder according to mode
		if mode == 'genome':
			BUSCO_report = HCGB_files.create_subfolder("BUSCO_assembly", outdir_report)
		elif mode == 'proteins':
			BUSCO_report = HCGB_files.create_subfolder("BUSCO_annot", outdir_report)

		## generate plots
		print ("+ Generate summarizing plots...")
		BUSCO_plots(dataFrame_results, BUSCO_report, options.threads)	
		print ('\n+ Check quality plots in folder: %s' %BUSCO_report)

		##	TODO 
		##	Parse BUSCO statistics in dataframe (stats_results) for discarding samples if necessary
		##	given a cutoff, discard or advise to discard some samples

		### print statistics
		stats_results.to_csv(BUSCO_report + "/BUSCO_stats.csv")
		name_excel = BUSCO_report + "/BUSCO_stats.xlsx"
		writer = pd.ExcelWriter(name_excel, engine='xlsxwriter')
		stats_results.to_excel(writer, sheet_name="BUSCO statistics")	
		writer.save()
		
		print ('\n+ Check quality statistics in folder: %s' %BUSCO_report)
	
	return(dataFrame_results)

################################################
def BUSCO_call(datasets, pd_samples, database_folder, threads, mode):
	## 
	## argument mode = proteins or genome
	##
	
	## get datasets
	print ("+ Check folder provided as database for available BUSCO datasets...")
	BUSCO_datasets = BUSCO_caller.BUSCO_retrieve_sets(datasets, database_folder)
	
	## BUSCO needs to chdir to output folder
	path_here = os.getcwd()
	
	print ("+ Checking quality for each sample retrieved...")
	## optimize threads: No need to optimize. There is a problem with the working dir of BUSCO and we need to change every time
	# We can use a with statement to ensure threads are cleaned up promptly
	with concurrent.futures.ThreadPoolExecutor(max_workers=1) as executor: ## need to do 1 by one as there is a problem with the working directory
		for DataSet in BUSCO_datasets:
			## send for each sample
			commandsSent = { executor.submit( BUSCO_runner, row['name'], DataSet, row['sample'], BUSCO_datasets[DataSet], row['busco_folder'], threads, mode): name for name, row in pd_samples.iterrows() }

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
	short_summary = pd.DataFrame(columns=('sample', 'dirname', 'name', 'ext', 'tag', 'busco_folder', 'busco_dataset', 'busco_summary', 'busco_results'))
	stats_summary = pd.DataFrame()

	## generate results
	for DataSet in BUSCO_datasets:
		for index, row in pd_samples.iterrows():
			#my_BUSCO_results_folder = row['busco_folder'] + '/run_' + DataSet
			my_BUSCO_results_folder = row['busco_folder'] + '/' + row['name'] + '/run_' + DataSet
			my_short_txt = 	my_BUSCO_results_folder + '/short_summary_' + DataSet + '.txt'
			
			if os.path.isfile(my_short_txt):
				short_summary.loc[len(short_summary)] = [ row['sample'], row['dirname'], row['name'], row['ext'], row['tag'], row['busco_folder'], DataSet, my_short_txt, my_BUSCO_results_folder ]
				my_stats = BUSCO_caller.BUSCO_stats(my_short_txt, row['name'], DataSet)
				stats_summary = pd.concat([stats_summary, my_stats])
		
	return (short_summary, stats_summary)

################################################
def BUSCO_runner(sample_name, DataSet, sample, dataset_path, output, threads, mode):
	#file_name = os.path.basename(sample)
	sample_name_dir = HCGB_files.create_subfolder(sample_name, output) ## to check in detached mode

	## run busco	
	code = BUSCO_caller.BUSCO_run( dataset_path, sample, threads, sample_name_dir, DataSet, mode)
	
	## retry it just in case
	if (code == 'FAIL'):
		BUSCO_caller.BUSCO_run( dataset_path, sample, threads, sample_name_dir, DataSet, mode)
	else:
		return ()

################################################
def BUSCO_plots(dataFrame_results, outdir, threads):

	## DataFrame columns ('sample', 'dirname', 'name', 'ext', 'tag', 'busco_folder', 'busco_dataset', 'busco_summary', 'busco_results'))
	list_datasets = set(dataFrame_results['busco_dataset'].tolist())
	list_samples = set(dataFrame_results['name'].tolist())

	plot_folder = HCGB_files.create_subfolder('BUSCO_plots', outdir)
	outdir_busco_plot = []
	
	print ("+ Get results for all samples summarized by dataset:")
	for dataset in list_datasets:
		print ("\t+ Get results for: ", dataset)
		plot_folder_dataset = HCGB_files.create_subfolder(dataset, plot_folder)
		outdir_busco_plot.append(plot_folder_dataset)
	
		for index, row in dataFrame_results.iterrows():
			if (dataset == row['busco_dataset']):
				shutil.copy(row['busco_summary'], plot_folder_dataset + '/short_summary_' + dataset + '_' + row['name'] + '.txt')
		
	print ("+ Get results for summarized by sample:")
	for sample in list_samples:
		print ("\t+ Get results for: ", sample)
		plot_folder_sample = HCGB_files.create_subfolder(sample, plot_folder)
		outdir_busco_plot.append(plot_folder_sample)

		for index, row in dataFrame_results.iterrows():
			if (sample == row['name']):
				shutil.copy(row['busco_summary'], plot_folder_sample + '/short_summary_' + row['busco_dataset'] + '_' + sample + '.txt')
	
	print ("+ Generate plots for each subset")
	
	## optimize threads
	threads_job = 1  ## threads optimization
	max_workers_int = threads

	# We can use a with statement to ensure threads are cleaned up promptly
	with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers_int) as executor: ## need to do 1 by one as there is a problem with the working directory
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

