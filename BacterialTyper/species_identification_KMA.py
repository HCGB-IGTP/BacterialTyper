#usr/bin/env python
'''
This code calls KMA software to find the best match (species identification) to the reads in one
or more fastq files or one fasta file in a (kmer) database produced using the KMA program (Philip T.L.C. Clausen, 
Frank M. Aarestrup & Ole Lund, "Rapid and precise alignment of raw reads against redundant databases with KMA", 
BMC Bioinformatics, 2018;19:307.) 
Jose F. Sanchez
Copyright (C) 2019 Lauro Sumoy Lab, IGTP, Spain
'''
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
from BacterialTyper import functions
from BacterialTyper import config
from BacterialTyper import ariba_caller

####################
#### DATABASE	####
####################

##################################################
def download_kma_database(folder, database, debug):

	## types: bacteria, archaea, protozoa, fungi, plasmids, typestrains
	## Default: downloads all "bacterial,plasmids" genomes from KMA website

	## ToDo: update with latest version, not working so far.
	#ftp_site = "ftp://ftp.cbs.dtu.dk/public/CGE/databases/KmerFinder/version/latest/"
	
	## ToDo: Set automatic: download config file and look for prefix for each sample and generate a dictionary to code the prefix for each db.
	
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
	
	## debug message
	if (debug):
		print (colored("Function call: download_kma_database " + folder + ' ' + database + '\n','yellow'))

	## prefix
	if (database == 'plasmids'):
		prefix = '.T'
	else:
		prefix = '.ATG'
		
	index_name = folder + '/' + database + prefix	

	## check if already download
	return_code_down = False
	if os.path.exists(folder):
		return_code_down = check_db_indexed(index_name)
		## debug message
		if (debug):
			print (colored("Folder database is already available:" + folder,'yellow'))
		
	if (return_code_down == False): ## folder does not exists

		## Download data
		print ("\t+ Downloading data now, it migth take a while....")

		ftp_site = "ftp://ftp.cbs.dtu.dk/public/CGE/databases/KmerFinder/version/20190107/"

		## debug message
		if (debug):
			print (colored("Download files via function wget_download:",'yellow'))

		url = ftp_site + database + '.tar.gz'
		#functions.wget_download(url, folder)

		md5_url = ftp_site + database + '.md5'
		#functions.wget_download(md5_url, folder)
		print ("\n\t+ Data succesfully downloaded.....")

		## get files
		files = os.listdir(folder)
		md5_sum = ""
		for f in files:
			if f.endswith('tar.gz'):
				tar_file = folder + '/' + f
			elif f.endswith('md5'):
				md5_sum = folder + '/' + f
		
		## check md5sum
		print ("\t+ Checking for integrity using md5sum")
		
		# get md5 sum from source
		md5_string = ""
		with open(md5_sum, 'r') as myfile:
			line = myfile.read()
		
		line = re.sub(r"\s", ',', line)
		md5_string = line.split(",")[0]
		
		## calculate md5 for file
		result_md5 = functions.check_md5sum(md5_string, tar_file)
		if (result_md5 == True):
		
			## debug message
			if (debug):
				print (colored("result md5sum matches code provided for file " + tar_file,'yellow'))

			# extract
			print ("\t+ Extracting database into destination folder: " + folder)
			functions.extract(tar_file, folder)	

		else:
			print (colored("*** ERROR: Some error ocurred during the downloading and file is corrupted ***", 'red'))
			return ("Error")
		
		## database should be unzipped and containing files...
		return_code_extract = check_db_indexed(index_name)
		
		if (return_code_extract):
			print("+ Database (%s) succesfully extracted in folder: %s..." %(database, folder))
		else:
			string = "*** ERROR: Some error ocurred during the extraction of the database (%s). Please check folder (%s) and downloading and file is corrupted ***" %(database, folder)
			print (colored(string, 'red'))
			return ("Error")

##################################################
def check_db_indexed(index_name):
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
			
	## dump in screen
	names = index_name + '.name'
	count = functions.get_number_lines(names)
	
	print ("\n\t+ Database seems OK and contains several entries (%s):\n" %count)

	if (count > 50):
		print ("\tToo many entries in the database.\n\tCheck file %s for further details." %names)
	else:
		names_hd = open(names, 'r')
		names_hd_read = names_hd.read()
		names_hd.close()
		print (names_hd_read)
	
	return(True)
	
##################################################
def index_database(fileToIndex, kma_bin, index_name, option):
	
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
	
	logFile = index_name + '.log'
	
	if (option == "new"):
		print ("\n+ Generate and index database for kmer alignment search...\n")
		cmd_kma_index = "%s index -i %s -o %s 2> %s" %(kma_bin, fileToIndex, index_name, logFile)
	elif (option == "add"):
		print ("\n+ Updating database with new entries...\n")
		cmd_kma_index = "%s index -i %s -o %s -t_db 2> %s" %(kma_bin, fileToIndex, index_name, logFile)

	functions.system_call(cmd_kma_index)	
	return_code = check_db_indexed(index_name)
	return(return_code)

########################
#### IDENTIFICATION	####
########################

##################################################
def kma_ident_call(out_file, files, sample_name, index_name, kma_bin, threads):
	###
	out_file_log = out_file + '.log'
	if len(files) == 2:
		#print ("Paired-end mode KMA search:\n")
		cmd_kma_search = "%s -Sparse -ipe %s %s -o %s -t_db %s -shm 1 -t %s 2> %s" %(kma_bin, files[0], files[1], out_file, index_name, threads, out_file_log)
	else:
		## to be tested
		print ("Single end mode KMA search:\n")
		cmd_kma_search = "%s -Sparse -i %s -o %s -t_db %s -shm 1 -t %s 2> %s" %(kma_bin, files[0], out_file, index_name, threads, out_file_log)

	return(functions.system_call(cmd_kma_search))
	
##################################################
def kma_ident_module(out_file, files, sample_name, index_name, threads):
	## kma_ident_call
	kma_bin = config.get_exe("kma")
	return(kma_ident_call(out_file, files, sample_name, index_name, kma_bin, threads))

##################################################
def parse_kma_results(out_file, cutoff):
	results = pd.read_csv(out_file, sep="\t")
	results_filter = results[results['Template_Coverage'] > cutoff] 
	## filter according to coverage of the template. 
	return (results_filter)

##################################################
def	help_options():
	print ("\nUSAGE: python %s file1 file2 name xxc threads path\n"  %os.path.realpath(__file__))

##################################################
def main():
	## this code runs when call as a single script

  	## control if options provided or help
	if len(sys.argv) > 1:
		print ("")
	else:
		help_options()
		exit()    	
	
	## to do: implement main function
	name = argv[1]
	cutoff = 80
	file_out = os.path.abspath(argv[2])
	parse_kma_results(file_out, cutoff)

##################################################
if __name__== "__main__":
	main()
	
	
	
	
################################################################################################
## 										KMA v1.2.2
################################################################################################	
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
	

