#!/usr/bin/env python3
##########################################################
## Jose F. Sanchez										##
## Copyright (C) 2019 Lauro Sumoy Lab, IGTP, Spain		##
##########################################################
'''
Calls Trimmomatic for the trimming of sequence adapter within fastq reads.
'''

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
from BacterialTyper import functions
from BacterialTyper import config

################################################
def	help_options():
	print ("\nUSAGE:\npython %s Outfolder file_R1 file_R2 trimmomatic threads trimmomatic_adapters sample_name\n"  %os.path.abspath(argv[0]))
	print ("\n*** If not paired-end provide 'na' for file_R2")

################################################
def trimmo_module(files, path_name, sample_name, threads, Debug, trimmomatic_adapters):
	## 
	## This functions generates a trimmomatic call using java and trimmomatic from 
	## the system with a minimun version (specified in config.py)
	## Checks if adapter file exists
	## Returns code from trimmo_call: OK/FAIL
	##
	
	## get exe
	trimmomatic_jar = config.get_exe('trimmomatic')
	java_path = config.get_exe('java')

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
	return(trimmo_call(java_path, path_name, sample_name, files, trimmomatic_jar, threads, trimmomatic_adapters, Debug))

################################################
def print_help_adapters():
	## [TODO]
	print (colored("\n\n***** TODO: Generate this help message *****\n\n", 'red'))

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

################################################
def trimmo_call(java_path, sample_folder, sample_name, files, trimmomatic_jar, threads, trimmomatic_adapters, Debug):
	##
	## Function to call trimmomatic using java. Can take single-end and pair-end files
	## sample_folder must exists before calling this function. 
	## It can be call from main or a module.
	## Returns code OK/FAIL according if suceeded or failed the system call
	## 

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

	## Paired or single end
	## set command
	cmd = ""	
	if (len(files) == 2): ## paired-end
		file_R1 = files[0]
		file_R2 = files[1]

		#print ('\t-', file_R2)
		trim_R1 = sample_folder + '/' + sample_name + '_trim_R1.fastq'
		orphan_R1 = sample_folder + '/' + sample_name + '_orphan_R1.fastq'
		trim_R2 = sample_folder + '/' + sample_name + '_trim_R2.fastq'
		orphan_R2 = sample_folder + '/' + sample_name + '_orphan_R2.fastq'

		cmd = "%s -jar %s PE -threads %s -trimlog %s %s %s %s %s %s %s ILLUMINACLIP:%s:2:30:10 LEADING:11 TRAILING:11 SLIDINGWINDOW:4:20 MINLEN:24 2> %s" %(java_path, trimmomatic_jar, threads, log_file, file_R1, file_R2, trim_R1, orphan_R1, trim_R2, orphan_R2, trimmomatic_adapters, trimmo_log)

	else: ## single end
		file_R1 = files[0]
		trim_R1 = sample_folder + '/' + sample_name + '_trim.fastq'

		cmd = "%s -jar %s SE -threads %s -trimlog %s %s %s ILLUMINACLIP:%s:2:30:10 LEADING:11 TRAILING:11 SLIDINGWINDOW:4:20 MINLEN:24 2> %s" %(java_path, trimmomatic_jar, threads, log_file, file_R1, trim_R1, trimmomatic_adapters, trimmo_log)

	## system call & return
	code = functions.system_call(cmd)
	if code == 'OK':
		## success stamps
		filename_stamp = sample_folder + '/.success'
		stamp =	functions.print_time_stamp(filename_stamp)	
		return('OK')	
	else:
		return('FAIL')	

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

	## call
	trimmo_call(path_name, sample_name, file_R1, file_R2, trimmomatic_jar, threads, trimmomatic_adapters)
		
################################################

'''******************************************'''
if __name__== "__main__":
	main()
		

