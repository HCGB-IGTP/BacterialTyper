#usr/bin/env python
'''
This code calls Trimmomatic for the trimming of sequence adapter within fastq reads.
Jose F. Sanchez
Copyright (C) 2019 Lauro Sumoy Lab, IGTP, Spain
'''

## useful imports
import time
import io
import os
import re
import subprocess
import sys
from sys import argv

## import my modules
pythonDir = os.path.dirname(os.path.realpath(__file__))
sys.path.append(pythonDir)
import functions

## import configuration
configDir = os.path.dirname(os.path.realpath(__file__)) + '/../../config/'
sys.path.append(configDir)
import config

######
def	help_options():
	print ("\nUSAGE:\npython %s Outfolder file_R1 file_R2 trimmomatic threads trimmomatic_adapters sample_name\n"  %os.path.abspath(argv[0]))

######	
def trimmo_module(file_R1, file_R2, path_name, sample_name, threads):
	trimmomatic_jar = config.EXECUTABLES('trimmomatic')
	trimmomatic_adapters = config.DATA('trimmomatic_adapters')
	## call
	return(trimmo_call(path_name, sample_name, file_R1, file_R2, trimmomatic_jar, threads, trimmomatic_adapters))

######
def trimmo_call(path_name, sample_name, file_R1, file_R2, trimmomatic_jar, threads, trimmomatic_adapters):
	
	print ("+ Cutting adapters for sample: ", sample_name)
	print ('\t-', file_R1)
	print ('\t-', file_R2)
		
	sample_folder = path_name + '/' + sample_name
	trimmo_log = path_name + '/' + sample_name + '.log'

	## create folder
	functions.create_folder(path_name)
	functions.create_folder(sample_folder)

	log_file = sample_folder + '/' + sample_name + '.log'
	trim_R1 = sample_folder + '/' + sample_name + '_trim_R1.fastq'
	trim_R2 = sample_folder + '/' + sample_name + '_trim_R2.fastq'
	orphan_R1 = sample_folder + '/' + sample_name + '_orphan_R1.fastq'
	orphan_R2 = sample_folder + '/' + sample_name + '_orphan_R2.fastq'
		
	cmd = "java -jar %s PE -threads %s -phred33 -trimlog %s %s %s %s %s %s %s ILLUMINACLIP:%s:2:30:10 LEADING:11 TRAILING:11 SLIDINGWINDOW:4:20 MINLEN:24 > %s" %(trimmomatic_jar, threads, log_file, file_R1, file_R2, trim_R1, orphan_R1, trim_R2, orphan_R2, trimmomatic_adapters, trimmo_log)
	return(functions.system_call(cmd))	

######
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
		
######

'''******************************************'''
if __name__== "__main__":
	main()
		

