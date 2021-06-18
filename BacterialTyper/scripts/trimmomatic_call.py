#!/usr/bin/env python3
##########################################################
## Jose F. Sanchez										##
## Copyright (C) 2019 Lauro Sumoy Lab, IGTP, Spain		##
##########################################################
"""
Calls Trimmomatic for the trimming of sequence adapter within fastq reads.
"""

## useful imports
import time
import io
import os
import re
import subprocess
import sys
from sys import argv
from termcolor import colored

## import my modules
import HCGB
import HCGB.functions.time_functions as HCGB_time
import HCGB.functions.system_call_functions as HCGB_sys

from BacterialTyper.config import set_config

################################################
def trimmo_module(files, path_name, sample_name, threads, Debug, trimmomatic_adapters, trimmomatic_params):
	## 
	## This functions generates a trimmomatic call using java and trimmomatic from 
	## the system with a minimum version (specified in config.py)
	## Checks if adapter file exists
	## Returns code from trimmo_call: OK/FAIL
	##
	
	## get exe
	trimmomatic_jar = set_config.get_exe('trimmomatic')
	java_path = set_config.get_exe('java')

	## check if it exists
	if os.path.isfile(trimmomatic_adapters):
		## debug message
		if (Debug):
			print (colored("**DEBUG: trimmomatic_adapters file exists **", 'yellow'))
			print (trimmomatic_adapters)
	else:
		## rise error & exit
		print (colored("***ERROR: Trimmomatic adapters file does not exist: " + trimmomatic_adapters,'red'))
		exit()
	
	## call
	return(trimmo_call(java_path, path_name, sample_name, files, trimmomatic_jar, threads, trimmomatic_adapters, trimmomatic_params, Debug))

################################################
def print_help_adapters():
	"""Trimmomatic adapters information
	
	BacterialTyper includes a file :file:`BacterialTyper.data.available_Trimmomatic_adapters.fasta` with default sequencing adapters provided by Trimmomatic_ (v0.39).
	
	User can provide using the option --adapters different sequencing adapters.
	
	.. seealso:: Additional information on Trimmomatic adapters available.
	
		- :doc:`Trimmomatic adapters <../../../data/trimmomatic_adapters>` 
	
	.. include:: ../../links.inc 
	"""
	## [TODO]
	print (colored("\n\n ***** TODO: Generate this help message *****\n\n", 'red'))
	
	print("BacterialTyper includes a file with default sequencing adapters provided by Trimmomatic (v0.39).")
	print("See additional details in the documentation.\n")
	
	print("Users can provide any sequencing adapter of interest using the option --adapters in fasta format file.")

################################################
def trimmo_call(java_path, sample_folder, sample_name, files, trimmomatic_jar, threads, trimmomatic_adapters, trimmomatic_params, Debug):
	##
	## Function to call trimmomatic using java. Can take single-end and pair-end files
	## sample_folder must exists before calling this function. 
	## It can be call from main or a module.
	## Returns code OK/FAIL according if succeeded or failed the system call
	## 

	#######################################
	## http://www.usadellab.org/cms/?page=trimmomatic
	#
	# ILLUMINACLIP:fasta_file.fa:2:30:10 LEADING:11 TRAILING:11 SLIDINGWINDOW:4:20 MINLEN:24
	#
	# This will perform the following:
	#	Remove adapters (ILLUMINACLIP:fasta_file.fa:2:30:10)
	#	Remove leading low quality or N bases (below quality 11) (LEADING:11)
	#	Remove trailing low quality or N bases (below quality 11) (TRAILING:11)
	#	Scan the read with a 4-base wide sliding window, cutting when the average quality per base drops below 20 (SLIDINGWINDOW:4:20)
	#	Drop reads below the 24 bases long (MINLEN:24)
	#######################################


	## debug message
	if (Debug):
		print (colored("+ Cutting adapters for sample: " + sample_name, 'yellow'))
		
	## log files
	log_file = sample_folder + '/' + sample_name + '_call.log'
	trimmo_log = sample_folder + '/' + sample_name + '.log'
	
	## init
	file_R1 = ""
	file_R2 = ""
	trim_R1 = ""
	orphan_R1 = ""
	trim_R2 = ""
	orphan_R2 = ""

	## conda installation includes a wrapper and no java jar call is required
	if trimmomatic_jar.endswith('jar'):
		cmd = "%s -jar %s"  %(java_path, trimmomatic_jar)
	else:
		cmd = "%s"  %(trimmomatic_jar)

	## Paired or single end
	## set command
	if (len(files) == 2): ## paired-end
		file_R1 = files[0]
		file_R2 = files[1]

		#print ('\t-', file_R2)
		trim_R1 = sample_folder + '/' + sample_name + '_trim_R1.fastq'
		orphan_R1 = sample_folder + '/' + sample_name + '_orphan_R1.fastq'
		trim_R2 = sample_folder + '/' + sample_name + '_trim_R2.fastq'
		orphan_R2 = sample_folder + '/' + sample_name + '_orphan_R2.fastq'

		cmd = cmd + " PE -threads %s -trimlog %s %s %s %s %s %s %s " %(threads, log_file, 
																	file_R1, file_R2, trim_R1, 
																	orphan_R1, trim_R2, orphan_R2)
	else: ## single end
		file_R1 = files[0]
		trim_R1 = sample_folder + '/' + sample_name + '_trim.fastq'

		cmd = cmd + " SE -threads %s -trimlog %s %s %s " %(threads, log_file, file_R1, trim_R1)

	## common parameters
	cmd = cmd + "ILLUMINACLIP:%s:%s LEADING:%s TRAILING:%s SLIDINGWINDOW:%s MINLEN:%s 2> %s" %(trimmomatic_adapters, 
																								trimmomatic_params['ILLUMINACLIP'],
																								trimmomatic_params['LEADING'],
																								trimmomatic_params['TRAILING'],
																								trimmomatic_params['SLIDINGWINDOW'],
																								trimmomatic_params['MINLEN'],
																								trimmo_log)
	## system call & return
	code = HCGB_sys.system_call(cmd)
	if code == 'OK':
		## success stamps
		filename_stamp = sample_folder + '/.success'
		stamp =	HCGB_time.print_time_stamp(filename_stamp)	
		return('OK')	
	else:
		return('FAIL')	


################################################
def	help_options():
	print ("\nUSAGE:\npython %s Outfolder file_R1 file_R2 trimmomatic threads trimmomatic_adapters sample_name\n"  %os.path.abspath(argv[0]))
	print ("\n*** If not paired-end provide 'na' for file_R2")

################################################
def main():

  	## control if options provided or help
	if len(sys.argv) > 1:
		print ("")
	else:
		help_options()
		exit()    	

	## ARGV
	path_name = os.path.abspath(argv[1])
	file_R1 = os.path.abspath(argv[2])
	file_R2 = os.path.abspath(argv[3])
	trimmomatic_jar = argv[4]
	threads = argv[5]
	trimmomatic_adapters = os.path.abspath(argv[6]) ##"/imppc/labs/lslab/share/data/references/Trimmomatic_adapters.fa"
	sample_name = argv[7]

	## default
	trimmomatic_params = {
		"ILLUMINACLIP": "2:30:10",
		"LEADING": "11",
		"TRAILING":"11",
		"SLIDINGWINDOW": "4:20",
		"MINLEN": "24"
		
		}

	## call
	trimmo_call(path_name, sample_name, file_R1, file_R2, 
			trimmomatic_jar, threads, trimmomatic_adapters, 
			trimmomatic_params, True)
		
################################################

'''******************************************'''
if __name__== "__main__":
	main()
		

