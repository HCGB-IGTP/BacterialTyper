#!/usr/bin/env python3
##########################################################
## Jose F. Sanchez										##
## Copyright (C) 2019 Lauro Sumoy Lab, IGTP, Spain		##
##########################################################
"""
Call KMA_ software to find the best match in reads file or fasta file in a (kmer) database produced using the KMA program 

.. include:: ../../links.inc 
"""
## useful imports
import time
import io
import os
import re
import sys
import pandas as pd
from sys import argv
from io import open
from termcolor import colored
import shutil

## import my modules
from BacterialTyper.config import set_config
from BacterialTyper.scripts import ariba_caller

import HCGB
import HCGB.functions.time_functions as HCGB_time
import HCGB.functions.system_call_functions as HCGB_sys
import HCGB.functions.files_functions as HCGB_files
import HCGB.functions.main_functions as HCGB_main

########################
#### INFORMATION	####
########################

##################################################
def help_kma_database():
	"""
	KMA_ software description message.
	"""
	## [TODO]
	print (colored("\n\n***** TODO: Generate this help message *****\n\n", 'red'))

####################
#### DATABASE	####
####################

##################################################
def download_kma_database(folder, database, debug):
	"""
	Downloads databases from KMA website.
	
	Using the latest available ftp datasets, this function downloads available datasets using
	function :func:`BacterialTyper.scripts.functions.wget_download`. 
	
	Ftp site: "ftp://ftp.cbs.dtu.dk/public/CGE/databases/KmerFinder/version/latest/"
	
	It also downloads the md5sum for the dataset selected and compares with the 
	
	:param folder: Absolute path to folder that contains database.
	:param database: Possible options: [bacteria, archaea, protozoa, fungi, plasmids, typestrains, viral].
	:param debug: True/false for printing debugging messages.
	
	:type folder: string
	:type database: string
	:type debug: boolean
	
	.. seealso:: This function depends on other ``BacterialTyper`` functions called:
	
		- :func:`BacterialTyper.scripts.functions.wget_download`

		- :func:`BacterialTyper.scripts.functions.check_md5sum`

		- :func:`BacterialTyper.scripts.functions.extract`
		
		- :func:`BacterialTyper.scripts.functions.print_time_stamp`
		
		- :func:`BacterialTyper.scripts.functions.read_time_stamp`

		- :func:`BacterialTyper.scripts.species_identification_KMA.check_db_indexed`

	"""

	## ToDo: update with latest version
	ftp_site = "http://www.cbs.dtu.dk/public/CGE/databases/KmerFinder/version/latest/"
	db_url = ftp_site + database + '.tar.gz'
	md5_url = ftp_site + database + '.md5'

	## In v20190107 there was a plasmid database.
	#ftp_site = "ftp://ftp.cbs.dtu.dk/public/CGE/databases/KmerFinder/version/20190107/"

	############################################################################
	## ToDo: Set automatic: download config file and look for prefix for each 
	## sample and generate a dictionary to code the prefix for each db.
	############################################################################
	
	# Database configuration file - Describes the content of the database
	# Each db consist of 5 files with the following extensions: b, comp.b, length.b, seq.b, name
	# Other important files are: .name, .kma.entries.all, .kma.entries.deleted, .kma.entries.added, .md5
	# db_prefix	name	description
	#bacteria.ATG	Bacteria Organisms	Bacteria organisms library prefix=ATG
	#plasmids.T	Bacteria Plasmids	Bacteria plasmids library prefix=T
	#typestrains.ATG	Bacteria Type Strains	Bacteria type strains library prefix=ATG
	#fungi.ATG	Fungi	Fungi library prefix=ATG
	#protozoa.ATG	Protozoa	Protozoa library prefix=ATG
	#archaea.ATG	Archaea	Archaea library prefix=ATG	
	
	HCGB_files.create_folder(folder)
	db_folder = HCGB_files.create_subfolder(database, folder)
	
	## debug message
	if (debug):
		print (colored("Function call: download_kma_database " + db_folder + ' ' + database + '\n','yellow'))

	## prefix
	if (database == 'plasmids'):
		prefix = '.T'
	elif (database == 'viral'):
		prefix = '.TG'
	else:
		prefix = '.ATG'
		
	index_name = os.path.join(db_folder, database + prefix)

	## check if already download
	return_code_down = False
	if os.path.exists(db_folder):
		return_code_down = check_db_indexed(index_name, db_folder)
		stamp =	HCGB_time.read_time_stamp(db_folder + '/.success')
		
		## debug message
		if (debug):
			print (colored("Folder database is already available:" + db_folder,'yellow'))
			
	if (return_code_down == False): ## folder does not exists

		## Download data
		print ("\t+ Downloading data now, it may take a while....")

		## debug message
		if (debug):
			print (colored("Download files via function wget_download:",'yellow'))
		
		## connect to db_url
		HCGB_sys.wget_download(db_url, db_folder)
		HCGB_sys.wget_download(md5_url, db_folder)
		print ("\n\t+ Data downloaded.....")

		## get files
		files = os.listdir(db_folder)
		md5_sum = ""
		for f in files:
			if f.endswith('tar.gz'):
				tar_file = db_folder + '/' + f
			elif f.endswith('md5'):
				md5_sum = db_folder + '/' + f
		
		## check md5sum
		print ("\t+ Checking for integrity using md5sum")
		
		# get md5 sum from source
		md5_string = ""
		with open(md5_sum, 'r') as myfile:
			line = myfile.read()
		
		line = re.sub(r"\s", ',', line)
		md5_string = line.split(",")[0]
		
		## calculate md5 for file
		result_md5 = HCGB_sys.check_md5sum(md5_string, tar_file) ## FIXME: Not conda supported
		if (result_md5 == True):
		
			## debug message
			if (debug):
				print (colored("result md5sum matches code provided for file " + tar_file,'yellow'))

			# extract
			print ("\t+ Extracting database into destination folder: " + db_folder)
			HCGB_files.extract(tar_file, db_folder)	

		else:
			print (colored("*** ERROR: Some error occurred during the downloading and file is corrupted ***", 'red'))
			return ("Error")
			
		## database should be unzipped and containing files...
		return_code_extract = check_db_indexed(index_name, db_folder)
		
		if (return_code_extract):
			print("+ Database (%s) successfully extracted in folder: %s..." %(database, db_folder))
		else:
			string = "*** ERROR: Some error occurred during the extraction of the database (%s). Please check folder (%s) and downloading and file is corrupted ***" %(database, folder)
			print (colored(string, 'red'))
			return ("Error")
		
		## print timestamp
		filename_stamp = db_folder + '/.success'
		stamp =	HCGB_time.print_time_stamp(filename_stamp)
		
	## create dictionary
	db_info = {'ftp_site':db_url,
			'database': database,
			'path': folder,
			'index_name': index_name,
			'time':stamp}
	return (db_info)


##################################################
def check_db_indexed(index_name, folder):
	"""
	Check the status of a database
	
	:param index_name: Index name for the database
	:param folder: Absolute path of the folder containing the database.
	
	:type index_name: string
	:type folder: string
	
	:returns: True/False for the index status.
	
	.. seealso:: This function depends on other ``BacterialTyper`` functions called:
		
		- :func:`BacterialTyper.scripts.functions.readList_fromFile`
		
		- :func:`BacterialTyper.scripts.functions.get_number_lines`
		
		- :func:`BacterialTyper.scripts.functions.read_time_stamp`

		- :func:`BacterialTyper.scripts.functions.print_time_stamp`
	
	"""
	
	# Each db consist of 5 files with the following extensions: b, comp.b, length.b, seq.b, name
	my_index_list = [".comp.b", ".index.b", ".length.b", ".name", ".seq.b"]

	print ("\t+ Checking if database has been previously indexed...")
	for sufix in my_index_list:
		##print (sufix)
		
		my_file = index_name + sufix
		if os.path.isfile(my_file):
			print ("\t" + my_file + ' exists...')
		else:
			if (sufix == '.index.b'):
				continue
			else:
				return(False)
	
	## check if previously assembled and succeeded
	filename_stamp = folder + '/.success'
	if os.path.isfile(filename_stamp):
		stamp =	HCGB_time.read_time_stamp(filename_stamp)
		print (colored("\tDatabase was generated on: %s" %stamp, 'yellow'))

		## Check if necessary to download again after several months/days
		days_passed = HCGB_time.get_diff_time(filename_stamp)
		print ("\t\t** %s days ago" %days_passed)		
		## download again
		if (days_passed > 60): 
			print ("\t\t** Downloading information again just to be sure...")
			return(False)
	
	## dump in screen
	names = index_name + '.name'
	count = HCGB_main.get_number_lines(names)
	
	print ("\n\t+ Database seems OK and contains several entries (%s):\n" %count)
	if (count > 50):
		print ("\tToo many entries in the database.\n\tCheck file %s for further details." %names)
	else:
		entries = HCGB_main.readList_fromFile(names)
		print (*entries, sep='\n')

	return(True)
	
##################################################
def index_database(fileToIndex, kma_bin, index_name, option, folder, type_option):
	"""
	Calls KMA_ software to index fasta files into a database for later KMA identification.
	
	:param fileToIndex: Fasta file to include in the database.
	:param kma_bin: Absolute path to kma executable binary. 
	:param index_name: Name for the database.
	:param option: Option to create or update a database.
	:param folder: Absolute path to folder containing database.
	:param type_option: Option to index the database with batch: batch [Default: Off].
	
	:type fileToIndex: string
	:type kma_bin: string 
	:type index_name: string 
	:type option: string 
	:type folder: string 
	:type type_option: string 	 
	
	:returns: It returns message from :func:`BacterialTyper.scripts.species_identification_KMA.check_db_indexed`.
	
	.. seealso:: This function depends on other ``BacterialTyper`` functions called:
	
		- :func:`BacterialTyper.scripts.functions.create_folder`
		
		- :func:`BacterialTyper.scripts.functions.system_call`
		
		- :func:`BacterialTyper.scripts.species_identification_KMA.check_db_indexed`
	"""	
	
	########################################################################################
	## 								KMA_index-1.2.2							
	########################################################################################
	# kma_index creates the databases needed to run KMA, from a list of fasta files given.
	# Options are:		
	#				Desc:									Default:
	#
	#	-i			Input/query file name (STDIN: "--")		None
	#	-o			Output file								Input/template file
	#	-batch		Batch input file
	#	-deCon		File with contamination (STDIN: "--")	None/False
	#	-batchD		Batch decon file
	#	-t_db		Add to existing DB						None/False
	#	-k			Kmersize								16
	#	-k_t		Kmersize for template identification	16
	#	-k_i		Kmersize for indexing					16
	#	-ML			Minimum length of templates	kmersize 	(16)	
	#	-CS			Start Chain size						1 M
	#	-ME			Mega DB									False
	#	-NI			Do not dump *.index.b					False
	#	-Sparse		Make Sparse DB ('-' for no prefix)		None/False
	#	-ht			Homology template						1.0
	#	-hq			Homology query							1.0
	#	-and		Both homolgy thresholds
	#				has to be reached						or
	#	-v			Version
	#	-h			Shows this help message
	#######################################################################################
	
	## check if file exists
	if os.path.isfile(index_name):
		index_file_name = index_name
	else:
		index_file_name = folder + '/' + index_name
		
	logFile = index_file_name + '.log'
	
	## check if folder exists
	HCGB_files.create_folder(folder)

	## single file
	if (type_option == 'batch'):
		type_option = '-batch'
	else:
		type_option = '-i'
	
	## new or add to existing db
	if (option == "new"):
		print ("\n+ Generate and index database for kmer alignment search...\n")
		cmd_kma_index = "%s index %s %s -o %s 2> %s" %(kma_bin, type_option, fileToIndex, index_file_name, logFile)
	elif (option == "add"):
		print ("\n+ Updating database with new entries...\n")
		cmd_kma_index = "%s index %s %s -o %s -t_db %s 2> %s" %(kma_bin, type_option, fileToIndex, index_file_name, index_file_name, logFile)

	code = HCGB_sys.system_call(cmd_kma_index)	
	if code == 'FAIL':
		print (colored("Database generated an error during the index: %s" %index_name, 'red'))
		print (colored("EXIT", 'red'))
		exit()
		
	return_code = check_db_indexed(index_file_name, folder)
	return(return_code)

##################################################
def generate_db(file_abs_paths, name, fold_name, option, type_option, Debug, kma_bin):
	"""Generate a call to create or update index KMA databases for later kmer identification. 

	:param file_abs_paths: List of absolute paths fasta genome files to include in the database.
	:param name: Database name.
	:param fold_name: Directory path to store database generated.
	:param option: Generate a new database (option = 'new') or add to pre-existing database (option = 'add'). If database exists, automatically adds.
	:param type_option: Index genome fasta files one by one (option_type='single') or using a batch file containing multiple entries (option='batch').
	:param kma_bin:	Binary executable for KMA software 
	:param Debug: True/False for debugging messages.
	
	:type file_abs_paths: list
	:type name: string
	:type fold_name: string
	:type option: string
	:type type_option: string 
	:type kma_bin:
	:type Debug: bool
		
	:returns: Absolute path to database generated
	
	.. seealso:: This function depends on other ``BacterialTyper`` functions called:
	
		- :func:`BacterialTyper.scripts.functions.readList_fromFile`
		
		- :func:`BacterialTyper.scripts.functions.printList2file`
		
		- :func:`BacterialTyper.scripts.species_identification_KMA.check_db_indexed`

		- :func:`BacterialTyper.scripts.species_identification_KMA.index_database`
		
	"""

	print ('+ Updating the KMA database: ', name)			

	## check
	if len(file_abs_paths) > 1:
		## read db in fold_name and get index files
		info = fold_name + '/' + name + '.db'

		## 
		lineList = []
		toIndexList = []
		indexedList = []		

		###
		if os.path.exists(info):
			lineList = HCGB_main.readList_fromFile(info)
			option = 'add'

		for f in file_abs_paths:
			baseName = os.path.basename(f)
			
			## check if already index
			if baseName in lineList:
				print (colored('\t+ File %s is already available in database %s' %(baseName, name), 'green'))
				indexedList.append(f)
			else:
				toIndexList.append(f)		
		
		if toIndexList:
			## generate batch and call
			info2 = fold_name + '/.batch_entries.txt'
			HCGB_main.printList2file(info2, toIndexList)
			status = index_database(info2, kma_bin, name, option, fold_name, type_option)
			final_list = set(lineList + toIndexList + indexedList)
			final_list_name = [os.path.basename(f) for f in final_list]
			HCGB_main.printList2file(info, final_list_name)
			count_files = len(toIndexList)
			print ('+ %s samples have been added to the database' %count_files)
		else:
			print ('\n+ No new sequences were added to the database.')
			return (fold_name + '/' + name)			
		
	else:
		file_name = file_abs_paths[0]
		## check if previously indexed
		status = check_db_indexed(file_name, fold_name)
		if (status): #true
			## debug message
			if (Debug):
				print (colored("**DEBUG: Database (%s) is indexed" %file_name + " **", 'yellow'))
			return (file_name)
		else: #false
			## debug message
			if (Debug):
				print (colored("**DEBUG: Database (%s) is not indexed" %file_name + " **", 'yellow'))
			status = index_database(file_name, kma_bin, file_name, option, fold_name, type_option)
	
	## return
	if (status): #true
		return (file_name)
	else:
		return False

def load_db(kma_bin, db2use):
	"""
	This function loads the given database in memory.
	
	:param kma_bin: Absolute path to KMA binary.
	:param db2use: Database to load in memory
	
	:type kma_bin: string
	:type db2use: string 
	
	:returns: System call status. 
	"""
	cmd_load_db = "%s shm -t_db %s -shmLvl 1" %(kma_bin, db2use)
	return_code_load = HCGB_sys.system_call(cmd_load_db)
	
	return(return_code_load)

def remove_db(kma_bin, db2use):
	"""
	This function removes the given database from memory.
	
	:param kma_bin: Absolute path to KMA binary.
	:param db2use: Database to remove from memory
	
	:type kma_bin: string
	:type db2use: string 
	
	:returns: System call status. 
	"""

	cmd_rm_db = "%s shm -t_db %s -shmLvl 1 -destroy" %(kma_bin, db2use)
	return_code_rm = HCGB_sys.system_call(cmd_rm_db)
	
	return(return_code_rm)
			

########################
#### IDENTIFICATION	####
########################

##################################################
def kma_ident_call(out_file, files, sample_name, index_name, kma_bin, option, threads):
	"""Create kma system call for kmer identification. 
	
	Paired-end end or single end fastq files accepted. It generates a time stamp if succeeds.

	:param out_file: Absolute path and basename for the output files generated with results.
	:param files: List of absolute paths for fastq files to search againts the database.
	:param sample_name: Directory path to store database generated.
	:param index_name: Database name
	:param kma_bin: Binary executable for KMA software.
	:param option: Additional options to pass to the system call.
	:param threads: Number of CPUs to use. 

	:type out_file: string
	:type files: list
	:type sample_name: string
	:type index_name: string
	:type kma_bin: string
	:type option: string
	:type threads: integer	

	:returns: System call returned finish status.
	
	.. seealso:: This function depends on other ``BacterialTyper`` functions called:
	
		- :func:`BacterialTyper.scripts.functions.system_call`
	
		- :func:`BacterialTyper.scripts.functions.print_time_stamp`

	"""

	###
	out_file_log = out_file + '.log'
	if len(files) == 2:
		cmd_kma_search = "%s -ipe %s %s -o %s -t_db %s -t %s %s 2> %s" %(kma_bin, files[0], files[1], out_file, index_name, threads, option, out_file_log)
	else:
		## TODO: test Single End
		cmd_kma_search = "%s -i %s -o %s -t_db %s -t %s %s 2> %s" %(kma_bin, files[0], out_file, index_name, threads, option, out_file_log)

	code = HCGB_sys.system_call(cmd_kma_search)

	if (code == 'OK'):
		## success stamps
		basename_tag = os.path.basename(out_file)
		folder = os.path.dirname(out_file)
		filename_stamp = folder + '.success_' + basename_tag
		stamp =	HCGB_time.print_time_stamp(filename_stamp)
		return('OK')
	else:
		return('FAIL')

##################################################
def parse_kma_results(out_file, cutoff):
	"""Filters KMA results given a cutoff threshold.
	
	:param out_file: Results file originated from kma system call (kma_ident_call)
	:param cutoff: Template coverage cutoff.

	:type out_file: string
	:type cutoff: integer	

	:returns: Filtered dataframe of results.
	"""
	results = pd.read_csv(out_file, sep="\t")
	results_filter = results[results['Template_Coverage'] > cutoff] 
	## filter according to coverage of the template. 
	return (results_filter)

##################################################
def	help_options():
	print ("\nUSAGE: python %s name fasta folder_name sample_name threads files_csv_string\n"  %os.path.realpath(__file__))

##################################################
def main():
	## this code runs when call as a single script

  	## control if options provided or help
	if len(sys.argv) > 1:
		print ("")
	else:
		help_options()
		exit()    	
	
	## arguments
	name = argv[1]
	fasta = os.path.abspath(argv[2])
	folder = os.path.abspath(argv[3])
	sample_name = argv[4]
	threads = argv[5]

	files = []
	for i,e in enumerate(argv):
		if i > 5:
			files.append(os.path.abspath(argv[i]))

	## other
	cutoff=80
	kma_bin = set_config.get_exe("kma")
	out_file = sample_name + ".out_kma-search.txt"

	## check if database is indexed
	if not  check_db_indexed(name, folder):
		index_database(fasta, kma_bin, name, 'new', folder, '')
		print ("\n+ Database indexed")
	

	## search files
	kma_ident_call(out_file, files, sample_name, folder + '/' + name, kma_bin, '', threads)

	##
	#parse_kma_results(file_out, cutoff)

##################################################
if __name__== "__main__":
	main()

################################################################################################
## 										KMA v1.2.2
################################################################################################	
## read additional information in: https://bitbucket.org/genomicepidemiology/kma/src/master/KMAspecification.pdf
# KMA-1.2.2 mapps raw reads to a template database.
# Options are:		
##				Desc:													Default:	Requirements:
#
#	-o			Output file												None		REQUIRED
#	-t_db		Template DB												None		REQUIRED
#	-i			Input file name(s)										STDIN
#	-ipe		Input paired end file name(s)
#	-int		Input interleaved file name(s)
#	-k			Kmersize												DB defined
#	-e			evalue													0.05
#	-ConClave	ConClave version										1
#	-mem_mode	Use kmers to choose best template, and save memory		False
#	-ex_mode	Searh kmers exhaustively								False
#	-ef			Print additional features								False
#	-vcf		Make vcf file, 2 to apply FT							False/0
#	-deCon		Remove contamination									False
#	-dense		Do not allow insertions in assembly						False
#	-ref_fsa	Consensus sequnce will have "n" instead of gaps			False
#	-matrix		Print assembly matrix									False
#	-a			Print all best mappings									False
#	-mp			Minimum phred score										20
#	-5p			Cut a constant number of nucleotides from the 5 prime.	0
#	-Sparse		Only count kmers										False
#	-Mt1		Map only to "num" template.								0 / False
#	-ID			Minimum ID												1.0%
#	-ss			Sparse sorting (q,c,d)									q
#	-pm			Pairing method (p,u,f)									u
#	-fpm		Fine Pairing method (p,u,f)								u
#	-apm		Sets both pm and fpm									u
#	-shm		Use shared DB made by kma_shm							0 (lvl)
#	-1t1		Skip HMM												False
#	-ck			Count kmers instead of pseudo alignment					False
#	-ca			Make circular alignments								False
#	-boot		Bootstrap sequence										False
#	-bc			Base calls should be significantly overrepresented.		[True]
#	-bc90		Base calls should be both  significantly 
#				overrepresented, and have 90% agreement.				False
#	-bcNano		Call bases at suspicious deletions, made for nanopore.	False
#	-bcd		Minimum depth at base									1
#	-bcg		Maintain insignificant gaps
#	-and		Both mrs and p_value thresholds has 
#				to reached to in order to report a template hit
#	-mq			Minimum mapping quality									0
#	-mrs		Minimum alignment score, normalized to 
#				alignment length										0.50
#	-reward		Score for match											1
#	-penalty	Penalty for mismatch									-2
#	-gapopen	Penalty for gap opening									-3
#	-gapextend	Penalty for gap extension								-1
#	-per		Reward for pairing reads								7
#	-cge		Set CGE penalties and rewards							False
#	-t			Number of threads										1
#	-v			Version
#	-h			Shows this help message
#################################################################################################

## kma index -batch ../database.txt -k 30 -k_t 30 -k_i 30 -o index_db_k30
## kma -ipe WTCHG_370809_205154/WTCHG_370809_205154_trim_R1.fastq WTCHG_370809_205154/WTCHG_370809_205154_trim_R2.fastq -o kma_search3 -t_db index_kma/index_db_k30 -t 3 -apm u -1t1 -k 30
	

