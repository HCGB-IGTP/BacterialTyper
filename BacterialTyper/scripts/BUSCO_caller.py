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

## import my modules
from BacterialTyper.scripts import functions
from BacterialTyper.config import set_config
from BacterialTyper import data

## TODO: append busco build lib to path
##sys.path.append(path2busco/busco/build/lib)

## TODO: append BUSCO Augustus path
## export AUGUSTUS_CONFIG_PATH=BUSCO/augustus_config/


###############
def busco_datasets():
	"""BUSCO dataset information
	
	:return: Dataframe containing information for each dataset available in file BUSCO_dataset.csv under data directory.
	
	Dataframe contains several fields: "Dataset", "Taxonomic range" and "ftp_site". Dataframe is indexed by Dataset. See an example:
	
	+-----------------+----------------------------------------------------------+----------------------------------------------------------------+
	| Dataset         | Taxonomic range                                          | ftp_site                                                       |
	+=================+==========================================================+================================================================+
	| deltaepsilonsub | phylum Proteobacteria – Delta and Epsilon proteobacteria | http://busco.ezlab.org/v2/datasets/deltaepsilonsub_odb9.tar.gz |
	+-----------------+----------------------------------------------------------+----------------------------------------------------------------+
	| bacteroidetes   | phylum Bacteroidetes                                     | http://busco.ezlab.org/v2/datasets/bacteroidetes_odb9.tar.gz   |
	+-----------------+----------------------------------------------------------+----------------------------------------------------------------+
	| tenericutes     + phylum Tenericutes                                       | http://busco.ezlab.org/v2/datasets/tenericutes_odb9.tar.gz     |
	+-----------------+----------------------------------------------------------+----------------------------------------------------------------+
	
	.. seealso:: Additional information on BUSCO datasets available.
	
		- :doc:`BUSCO datasets <../../../data/BUSCO_datasets>` 
	
	.. seealso:: This function depends on other BacterialTyper functions called:
	
		- :func:`BacterialTyper.data.data_files.data_list`
		
		- :func:`BacterialTyper.scripts.functions.file2dataframe`

	
	"""
	## read from file: BUSCO_dataset.csv
	BUSCO_dataset_file = data.data_files.data_list("BUSCO_dataset")
	busco_data_columns = ["Taxonomic range","Dataset","ftp_site"]
	busco_data = functions.file2dataframe(BUSCO_dataset_file, names=busco_data_columns)
	busco_data = busco_data.set_index('Dataset')
	return(busco_data)

###############
def print_help_BUSCO():
	functions.print_sepLine("*", 50, 'yellow')
	print ('BUSCO Help')
	functions.print_sepLine("*", 50, 'yellow')
	print ("Benchmarking of Universal Single Copy Orhtologs (BUSCO)\n")	
	
	print ("BUSCO provides quantitative measures for the assessment of genome assembly, gene set, and transcriptome completeness.")
	print ("It is based on evolutionary-informed expectations of gene content from near-universal single-copy orthologs selected from OrthoDB v9.\n")
	
	print ("BUSCO assessments are implemented in open-source software, with a large selection of lineage-specific sets of Single-Copy Orthologs.")
	print ("These conserved orthologs are ideal candidates for large-scale phylogenomics studies, and the annotated BUSCO gene models built during ")
	print ("genome assessments provide a comprehensive gene predictor training set for use as part of genome annotation pipelines.\n")
	
	print ("For more information visit the website: https://busco.ezlab.org/\n")
	functions.print_sepLine("*", 50, False)

	print ("\nPlease select among the available datasets according to your sample expected taxonomic range.") 
	print ("Bear in mind that several datasets could be provided as some represent a broad taxonomic range: ")
	print ("E.g. Bacteria > Proteobacteria > Gammaproteobacteria > Enterobacteriales")
	print ("E.g. Bacteria > Firmicutes > Bacillales")
	print ("Datasets:\n")	
	print_available_BUSCO()

###############
def print_available_BUSCO():
	print_df = busco_datasets()
	print_df = print_df.set_index('Dataset')
	pd.set_option('display.max_colwidth', None)
	functions.print_sepLine("-", 100, False)
	print (print_df)
	functions.print_sepLine("-", 100, False)
	print ("\n")

###############
def BUSCO_download(name, ftp, folder):
	"""Downloads BUSCO datasets provided
	
	This code checks if dataset is already available in folder provided. If not available proceeds to download it.
	
	It creates a subfolder (using :func:`BacterialTyper.scripts.functions.create_folder`) and downloads the dataset from ftp site provided (using :func:`BacterialTyper.scripts.functions.wget_download`).

	The file donwloaded would gunzipped so it is decompressed (using  :func:`BacterialTyper.scripts.functions.extract`). 
	
	A timestamp is printed to reflect download data (using :func:`BacterialTyper.scripts.functions.print_time_stamp`).
	
	The data download is checked for integrity using the function :func:`BacterialTyper.scripts.BUSCO_caller.BUSCO_check_dataset`.

	If the process failed, the whole process is retried (CAUTION: it might generate infinite loop).	
	
	:param name: Dataset name provided.
	:param ftp: FTP site for dataset.
	:param folder: Absolute path to folder to store results (e.g. database/BUSCO).
	:type name: string
	:type ftp: string
	:type folder: string
	
	:returns: Folder absolute path containing downloaded dataset.
	:rtype: string
	:warnings: Returns **FAIL** if check process failed.
	
	
	.. seealso:: This function depends on other BacterialTyper functions called:
	
		- :func:`BacterialTyper.scripts.functions.create_folder`
		
		- :func:`BacterialTyper.scripts.functions.wget_download`
		
		- :func:`BacterialTyper.scripts.functions.extract`
		
		- :func:`BacterialTyper.scripts.functions.print_time_stamp`
		
		- :func:`BacterialTyper.scripts.functions.print_sepLine`

		- :func:`BacterialTyper.scripts.BUSCO_caller.BUSCO_check_dataset`

	"""
	
	print (colored("\n+ BUSCO dataset: " + name + " - v9 OrthoDB", 'yellow')) ## modify OrthoDB version if changed
	subfolder = folder + '/' + name
	file_name = os.path.basename(ftp)
	path_file = subfolder + '/' + file_name
	folderName = subfolder + '/' + file_name.split('.tar.gz')[0]

	if os.path.exists(subfolder):
		print ('Already available in the path provided...') 
	else:
		## does not exists
		functions.create_folder(subfolder)

		functions.wget_download(ftp, subfolder)

		## extract
		print ("+ Extract file...")
		functions.extract(path_file, subfolder)

		## timestamp
		filename_stamp = subfolder + '/.success'
		functions.print_time_stamp(filename_stamp)

	## check if it is allright
	code = BUSCO_check_dataset(folderName)
	functions.print_sepLine("-", 50, False)
	
	if (code == 'FAIL'):
		print (colored('*** Dataset failed. Try to download it again...','red'))
		shutil.rmtree(folder)
		BUSCO_download(name, ftp, folder)
		
		### CAUTION: check for infinite loop
		
	return (folderName)

###############
def BUSCO_check_dataset(folder):
	config_file = folder + '/dataset.cfg'
	##	name=bacteria_odb9
	##	species=E_coli_K12
	##	domain=prokaryota
	##	creation_date=2016-11-01
	##	number_of_BUSCOs=148
	##	number_of_species=3663

	#print ("+ Checking the integrity of BUSCO dataset in folder: ", folder)
	functions.print_sepLine("+", 10, False)
	print ("Statistics")
	functions.print_sepLine("+", 10, False)
	if os.path.isfile(config_file):
		list_config = functions.readList_fromFile(config_file)
		for elem in list_config:
			line = elem.split("=")
			line[0] = line[0].replace("_", " ")
			print ("\t".join(line))
		
		if os.path.isfile(folder + '/../.success'):
			timestamp = functions.read_time_stamp(folder + '/../.success')
			print ("download date: ", timestamp)
		
		print ("Available in folder: ", folder)
		print (colored("Dataset....[ OK ]\n", 'green'))	
	else:
		print (colored("Dataset....[ FAIL ]\n", 'red'))	
		return ('FAIL')

###############
def BUSCO_retrieve_sets(list_datasets, folder):
	"""Retrieves datasets information available
	
	This functions checks in 'folder' if datasets in 'list_datasets' are available or downloads them if necessary. 
	Retrieves information using function :func:`BacterialTyper.scripts.busco_caller.busco_datasets` and checks if it is available or downloads it using :func:`BacterialTyper.scripts.BUSCO_caller.BUSCO_download`. 

	:param list_datasets: list of datasets of interest to check.
	:param folder: absolute path to folder that will contain all datasets of interest.
	 
	:type list_datasets: list
	:type folder: string
	
	:returns: Dictionary containing for each dataset of interest (key), the absolute path to the folder containing information downloaded (value).
	"""

	## check if name matches
	data_df = busco_datasets()
	for elem in list_datasets:
		if not elem in data_df.index:
			print (colored("***ATTENTION: Please check the datasets provided: " + elem, 'red'))
			print ("\n\nAvailable datasets:")
			print_available_BUSCO()
			exit()
	
	## database: folder
	dataset = {}
	for elem in list_datasets:
		ftp = data_df.loc[elem]['ftp_site']
		dataset[elem] = BUSCO_download(elem, ftp, folder)
	
	return (dataset)

###############
def BUSCO_run(dataset, fasta, threads, output_name, dataset_name, mode):

	my_out_folder = output_name + '/run_' + dataset_name
	## timestamp
	filename_stamp =  my_out_folder + '/.success'

	print (colored("\tBUSCO Dataset [%s]; Sample [%s]" %(dataset_name, fasta), 'yellow'))
		
	## check previous run
	if os.path.isfile(filename_stamp):
		timestamp = functions.read_time_stamp(filename_stamp)
		print (colored("\tSuccessfully run on date: %s"  %timestamp, 'green'))
	else:
	
		busco_bin = set_config.get_exe('busco')
		os.chdir(output_name)
		logFile = dataset_name + '.log'
		cmd = '%s -i %s -f -c %s --blast_single_core --mode %s -l %s -o %s > %s' %(busco_bin, fasta, threads, mode, dataset, dataset_name, logFile)
		functions.system_call(cmd)
	
		if os.path.isfile(my_out_folder + '/short_summary_' + dataset_name + '.txt'):
			## timestamp
			functions.print_time_stamp(filename_stamp)
		else:
			print (colored("BUSCO failed: Dataset [%s]; Sample [%s]" %(dataset_name, fasta), 'red'))
			return ('FAIL')

	return()

###############
def BUSCO_stats(summary_file, sample, dataset):
	lines_file = functions.get_info_file(summary_file)
	list_entries_num = []
	for l in lines_file:
		## discard messages and blank lines
		if not l.startswith("#") and l.strip():		
			## print each line
			## print (l)

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
				name_search = re.search(r".*\((.*)\)$", entries2[2])
				if name_search:
					pct_group = list_entries_dict[ name_search.group(1) ]					

					number = entries2[1] + ' (' + str(pct_group) + '%)'
					list_entries_num.append((entries2[2],  number))
				else:
					list_entries_num.append((entries2[2], entries2[1]))
					
	
	list_entries_num = list_entries_num + list_entries_pct
	stats = pd.DataFrame(list_entries_num, columns=('Type', sample))

	stats_returned = stats.set_index('Type').transpose()
	stats_returned['Database'] = dataset
	return (stats_returned)	

###############
def BUSCO_plot(outfolder):
	busco_plot_bin = set_config.get_exe('busco_plot')
	#logFile = dataset_name + '.log'
	cmd = '%s -wd %s' %(busco_plot_bin, outfolder)
	functions.system_call(cmd)

###############
def main():
	## this code runs when call as a single script

  	## control if options provided or help
	if not len(sys.argv) > 1:
		print_help_BUSCO()
		exit()
	
	BUSCO_stats(sys.argv[1], 'example', 'dataset1')
	
######
'''******************************************'''
if __name__== "__main__":
	main()
		

