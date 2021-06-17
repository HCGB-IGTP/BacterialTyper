#!/usr/bin/env python3
##########################################################
## Jose F. Sanchez										##
## Copyright (C) 2019 Lauro Sumoy Lab, IGTP, Spain		##
##########################################################
"""
Calls BUSCO software for quality control of annotation and assembly datasets
"""
## useful imports
import time
import io
import os
import re
import sys
import pandas as pd
import shutil
from termcolor import colored
import concurrent.futures

## import my modules
import HCGB
import HCGB.functions.time_functions as HCGB_time
import HCGB.functions.system_call_functions as HCGB_sys
import HCGB.functions.aesthetics_functions as HCGB_aes
import HCGB.functions.files_functions as HCGB_files
import HCGB.functions.main_functions as HCGB_main

from BacterialTyper.config import set_config
from BacterialTyper import data

##############################
def print_help_BUSCO():
	HCGB_aes.print_sepLine("*", 50, 'yellow')
	print ('BUSCO Help')
	HCGB_aes.print_sepLine("*", 50, 'yellow')
	print ("Benchmarking of Universal Single Copy Orthologs (BUSCO)\n")	
	
	print ("BUSCO provides quantitative measures for the assessment of genome assembly, gene set, and transcriptome completeness.")
	print ("It is based on evolutionary-informed expectations of gene content from near-universal single-copy orthologs selected from OrthoDB v9.\n")
	
	print ("BUSCO assessments are implemented in open-source software, with a large selection of lineage-specific sets of Single-Copy Orthologs.")
	print ("These conserved orthologs are ideal candidates for large-scale phylogenomics studies, and the annotated BUSCO gene models built during ")
	print ("genome assessments provide a comprehensive gene predictor training set for use as part of genome annotation pipelines.\n")
	
	print ("For more information visit the website: https://busco.ezlab.org/\n")
	HCGB_aes.print_sepLine("*", 50, False)

	print ("\nPlease select among the available datasets according to your sample expected taxonomic range.") 
	print ("Bear in mind that several datasets could be provided as some represent a broad taxonomic range: ")
	print ("E.g. Bacteria > Proteobacteria > Gammaproteobacteria > Enterobacteriales")
	print ("E.g. Bacteria > Firmicutes > Bacillales")
	print ("Datasets:\n")	
	print_available_BUSCO()

##############################
def print_available_BUSCO():
	HCGB_aes.print_sepLine("-", 100, False)
	busco_bin = set_config.get_exe('busco')
	
	## get datasets
	busco_bin_call = busco_bin + ' --list-datasets > tmp'
	HCGB_sys.system_call(busco_bin_call, message=False)
	
	## dump in screen
	with open("./tmp", 'r') as f:
		print(f.read())	
	
	## clean
	list_files = HCGB_main.get_fullpath_list("./busco_downloads")
	list_files + ['tmp']
	for i in list_files:
		os.remove(i)
	os.rmdir("./busco_downloads/information")		
	os.rmdir("./busco_downloads/")
		
	HCGB_aes.print_sepLine("-", 100, False)
	print ("\n")

##############################
def busco_datasets():
	"""BUSCO dataset information
	
	:return: List containing information for each dataset available in file BUSCO_dataset.csv under data directory.
	
	.. seealso:: Additional information on BUSCO datasets available.
	
		- :doc:`BUSCO datasets <../../../data/BUSCO_datasets>` 
	
	.. seealso:: This function depends on other BacterialTyper functions called:
	
		- :func:`BacterialTyper.data.data_files.data_list`
		
		- :func:`BacterialTyper.scripts.functions.file2dataframe`
	"""
	## read from file: BUSCO_dataset.csv
	BUSCO_dataset_file = data.data_files.data_list("BUSCO_dataset")
	busco_data = HCGB_main.get_data(BUSCO_dataset_file, ",", options="")
	return(busco_data)

##############################
def BUSCO_retrieve_sets(list_datasets, folder):
	"""Retrieves datasets information available
	
	This functions checks if list of datasets are correctle named and if there are already available in 'folder'. 
	Retrieves information using function :func:`BacterialTyper.scripts.busco_caller.busco_datasets` 
	and checks if it is already available and correct. 

	:param list_datasets: list of datasets of interest to check.
	:param folder: absolute path to folder that will contain all datasets of interest.
	 
	:type list_datasets: list
	:type folder: string
	
	:returns: Dictionary containing for each dataset of interest (key), the absolute path to the folder containing information downloaded (value).
	"""

	## check auto-lineage is alone
	if "auto-lineage" in list_datasets:
		if len(list_datasets) > 1:
			print ("\nDatasets provided:")
			print("\t".join(list_datasets))
			print (colored("***ATTENTION: It is not possible to provide Auto-lineage and a dataset at the same time: ", 'red'))
			exit()
	

	## check if name matches
	data_df = busco_datasets()
	for elem in list_datasets:
		exit_1 = True
		for column in data_df:
			if elem in data_df[column].values:
				exit_1=False
		## 
		if exit_1:
			if not elem == "auto-lineage":
				print (colored("***ATTENTION: \nNot available dataset. Please check the dataset provided: " + elem, 'red'))
				print ("\n\nAvailable datasets:")
				print_available_BUSCO()
				exit()
			else:
				print (colored("Auto-lineage mode provided:...... OK", 'green'))
		else:
			print (colored("Dataset provided: " + elem + "...... OK", 'green'))

	## check if in database: folder
	dataset = {}
	for elem in list_datasets:
		if not elem == "auto-lineage":
			## folder + lineages + elem
			status = BUSCO_check_dataset(os.path.abspath( os.path.join (folder, 'lineages', elem)), elem)
	
	return (list_datasets)

##############################
def BUSCO_check_dataset(folder, name):
	config_file = folder + '/dataset.cfg'
	##	name=bacteria_odb9
	##	species=E_coli_K12
	##	domain=prokaryota
	##	creation_date=2016-11-01
	##	number_of_BUSCOs=148
	##	number_of_species=3663

	if os.path.isdir(folder):
		#print ("+ Checking the integrity of BUSCO dataset in folder: ", folder)
		HCGB_aes.print_sepLine("+", 10, False)
		print ("Statistics for dataset: ")
		HCGB_aes.print_sepLine("+", 10, False)
		if os.path.isfile(config_file):
			list_config = HCGB_main.readList_fromFile(config_file)
			for elem in list_config:
				line = elem.split("=")
				line[0] = line[0].replace("_", " ")
				print (" "+ "\t".join(line))
			
			print()		
			print ("Available in folder: ", folder)
			print (colored("Dataset....[ OK ]\n", 'green'))	
		else:
			print (colored("Dataset....[ FAIL ]\n", 'red'))
			print ("+ Removing dataset to avoid further errors:")
			os.rmdir(folder)	
			return ('FAIL')

##############################
def BUSCO_runner(sample_name, sample, DataSet, output, threads, mode, database_folder):

	## run busco	
	code = BUSCO_run(sample_name, sample, threads, output, DataSet, mode, database_folder)
	
	## retry it just in case
	if (code == 'FAIL'):
		BUSCO_run(sample_name, sample, threads, output, DataSet, mode, database_folder)
	else:
		return ()
	
##############################
def BUSCO_call(datasets, pd_samples, database_folder, threads, mode):
	## 
	## argument mode = proteins or genome
	##
	
	## get datasets
	print ("+ Check folder provided as database for available BUSCO datasets...")
	BUSCO_datasets = BUSCO_retrieve_sets(datasets, database_folder)
	
	## ATTENTION: BUSCO needs to chdir to output folder
	path_here = os.getcwd()
	
	print ("+ Checking quality for each sample retrieved...")
	
	## optimize threads: No need to optimize. There is a problem with the working dir of BUSCO and we 
	## need to change every time. We can use a with statement to ensure threads are cleaned up promptly
	with concurrent.futures.ThreadPoolExecutor(max_workers=1) as executor: ## need to do 1 by one as there is a problem with the working directory
		for DataSet in BUSCO_datasets:
			## send for each sample
			commandsSent = { executor.submit( BUSCO_runner, row['name'],
											row['sample'], DataSet, 
											row['busco_folder'], threads, mode, database_folder): name for name, row in pd_samples.iterrows() }
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
	short_summary = pd.DataFrame(columns=('sample', 'dirname', 'name', 'ext', 'tag', 'busco_folder', 
										'busco_dataset', 'busco_summary', 'busco_results'))
	stats_summary = pd.DataFrame()

	## generate results
	for DataSet in BUSCO_datasets:
		for index, row in pd_samples.iterrows():
			#my_BUSCO_results_folder = row['busco_folder'] + '/run_' + DataSet
			my_BUSCO_results_folder = row['busco_folder'] + '/' + DataSet + '/run_' + DataSet
			my_short_txt = 	my_BUSCO_results_folder + '/short_summary.txt'
			
			if os.path.isfile(my_short_txt):
				short_summary.loc[len(short_summary)] = [ row['sample'], row['dirname'], row['name'], row['ext'], row['tag'], row['busco_folder'], DataSet, my_short_txt, my_BUSCO_results_folder ]
				my_stats = BUSCO_stats(my_short_txt, row['name'], DataSet)
				stats_summary = pd.concat([stats_summary, my_stats])
		
	#print(short_summary)
	#print(stats_summary)
	return (short_summary, stats_summary)

##############################
def BUSCO_run(sample_name, fasta, threads, output_name, dataset_name, mode, busco_db):

	my_out_folder = os.path.join(output_name, dataset_name + '/run_' + dataset_name)
	## timestamp
	filename_stamp =  my_out_folder + '/.success'

	print (colored("\tBUSCO Dataset [%s]; Sample [%s]" %(dataset_name, sample_name), 'yellow'))
		
	## check previous run
	if os.path.isfile(filename_stamp):
		timestamp = HCGB_time.read_time_stamp(filename_stamp)
		print (colored("\tSuccessfully run on date: %s"  %timestamp, 'green'))
	else:
	
		busco_bin = set_config.get_exe('busco')
		os.chdir(output_name)
		
		## init cmd configuration
		cmd = '%s -f -i %s -c %s --mode %s --download_path %s ' %(busco_bin, fasta, threads, mode, busco_db)
		
		## options if autolineage or given dataset
		if "auto-lineage" == dataset_name:
			logFile = 'auto_lineage.log'
			cmd = cmd + '--auto-lineage -o %s > %s' %(dataset_name, logFile)
		else:
			logFile = dataset_name + '.log'
			cmd = cmd + '-l %s -o %s > %s' %(dataset_name, dataset_name, logFile)
		
		## system call
		HCGB_sys.system_call(cmd)
		
		if os.path.isfile(my_out_folder + '/short_summary.txt'):
			## timestamp
			HCGB_time.print_time_stamp(filename_stamp)
		else:
			print (colored("BUSCO failed: Dataset [%s]; Sample [%s]" %(dataset_name, fasta), 'red'))
			return ('FAIL')

	return()

##############################
def BUSCO_stats(summary_file, sample, dataset):
	lines_file = HCGB_main.get_info_file(summary_file)
	list_entries_num = []
	for l in lines_file:
		## discard messages and blank lines
		if not l.startswith("#") and l.strip():		
			init_search = re.search(r"\*.*", l)
			version_search = re.search(r".*version.*", l)
			if not init_search and not version_search:
				#print each line
				#print (l)
				if l.startswith("\tC:"):
					l = l.replace('\t', '')
					l = l.replace('[', ',')
					l = l.replace(']', '')
					l = l.replace('%', '')
					entries = l.split(',')
					
					list_entries_pct = []
					list_entries_dict = {}
					for i in entries:
						list_entries_pct.append((i.split(":")))
						list_entries_dict[i.split(":")[0]] = i.split(":")[1]
		
				else:
					
					entries2 = l.split('\t')
					if not (entries2[1].startswith('hmm')) and not (entries2[1].startswith('prod')):
						name_search = re.search(r".*\((.*)\)$", entries2[2])
						
						if name_search:
							pct_group = list_entries_dict[ name_search.group(1) ]					
		
							number = entries2[1] + ' (' + str(pct_group) + '%)'
							list_entries_num.append((entries2[2],  number))
						else:
							list_entries_num.append((entries2[2], entries2[1]))
						
				
	#print(list_entries_dict)
	list_entries_num = list_entries_num + list_entries_pct
	stats = pd.DataFrame(list_entries_num, columns=('Type', sample))

	stats_returned = stats.set_index('Type').transpose()
	stats_returned['Database'] = dataset
	return (stats_returned)	

###############
def BUSCO_plot(outfolder):
	busco_plot_bin = set_config.get_exe('generate_plot')
	
	os.chdir(outfolder)
	#logFile = dataset_name + '.log'
	cmd = '%s -wd %s' %(busco_plot_bin, outfolder)
	HCGB_sys.system_call(cmd)
	return()
	
################################################

################################################
def BUSCO_plots(dataFrame_results, outdir, threads):

	## DataFrame columns ('sample', 'dirname', 'name', 'ext', 'tag', 'busco_folder', 'busco_dataset', 'busco_summary', 'busco_results'))
	list_datasets = set(dataFrame_results['busco_dataset'].tolist())
	list_samples = set(dataFrame_results['name'].tolist())

	plot_folder = HCGB_files.create_subfolder('BUSCO_plots', outdir)
	outdir_busco_plot = []
	
	## summary for dataset
	print ("+ Get results for all samples summarized by dataset:")
	for dataset in list_datasets:
		print ("\t+ Get results for: ", dataset)
		plot_folder_dataset = HCGB_files.create_subfolder(dataset, plot_folder)
		outdir_busco_plot.append(plot_folder_dataset)
	
		for index, row in dataFrame_results.iterrows():
			if (dataset == row['busco_dataset']):
				shutil.copy(row['busco_summary'], 
						plot_folder_dataset + '/short_summary.specific.' + dataset + '.' + row['name'] + '.txt')
		
	print ("+ Get results for summarized by sample:")
	for sample in list_samples:
		print ("\t+ Get results for: ", sample)
		plot_folder_sample = HCGB_files.create_subfolder(sample, plot_folder)
		outdir_busco_plot.append(plot_folder_sample)

		for index, row in dataFrame_results.iterrows():
			if (sample == row['name']):
				shutil.copy(row['busco_summary'], 
						plot_folder_sample + '/short_summary.specific.' + dataset + '.' + row['name'] + '.txt')
						#plot_folder_sample + '/short_summary.' + row['busco_dataset'] + '.' + sample + '.txt')
	
	print ("+ Generate plots for each subset")
	path_here = os.getcwd()
	
	for plot in outdir_busco_plot:
		BUSCO_plot(plot) 
			
	print ("+ All plots generated...")
	print ("+ Check results under folders in : ", plot_folder)
	
	os.chdir(path_here)
	return()
		

###############
def main():
	## this code runs when call as a single script

  	## control if options provided or help
	if not len(sys.argv) > 1:
		print_help_BUSCO()
		exit()
	
	list_datasets = str(sys.argv[1]).split(',')
	print (list_datasets)
	list_dataset =  BUSCO_retrieve_sets(list_datasets, sys.argv[2])
	
	for l in list_dataset:
		outp = os.path.join(sys.argv[4], l)
		BUSCO_run(os.path.abspath(sys.argv[3]), 1, sys.argv[4], l, "genome", os.path.abspath(sys.argv[2]))
	
######
'''******************************************'''
if __name__== "__main__":
	main()
		