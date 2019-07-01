 #!/usr/bin/env python3
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
from BacterialTyper import functions
from BacterialTyper import config

######
def	help_options():
	print ("\nUSAGE:\npython %s Outfolder file_R1 file_R2 trimmomatic threads trimmomatic_adapters sample_name\n"  %os.path.abspath(argv[0]))
	print ("\n*** If not paired-end provide 'na' for file_R2")

######	
def trimmo_module(files, path_name, sample_name, threads, Debug, trimmomatic_adapters):
	
	#trimmomatic_jar = "/software/debian-8/bio/trimmomatic-0.36/trimmomatic.jar" #config.EXECUTABLES['trimmomatic'] ## check if it works
	trimmomatic_jar = config.get_exe('trimmomatic')
	
	## check if it exists
	#trimmomatic_adapters = config.DATA['trimmomatic_adapters']
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
	return(trimmo_call(path_name, sample_name, files, trimmomatic_jar, threads, trimmomatic_adapters, Debug))

######
def trimmo_call(path_name, sample_name, files, trimmomatic_jar, threads, trimmomatic_adapters, Debug):
	
	##
	java_path = config.get_exe('java')

	## debug message
	if (Debug):
		print (colored("+ Cutting adapters for sample: " + sample_name, 'yellow'))
		
	## create folder
	sample_folder = path_name + '/' + sample_name
	functions.create_folder(path_name)
	functions.create_folder(sample_folder)

	log_file = sample_folder + '/' + sample_name + '.log'
	trimmo_log = path_name + '/' + sample_name + '.log'
	
	file_R1 = ""
	file_R2 = ""
	trim_R1 = ""
	orphan_R1 = ""
	trim_R2 = ""
	orphan_R2 = ""

	## Paired or single end
	if (len(files) == 2):
		file_R1 = files[0]
		file_R2 = files[1]

		#print ('\t-', file_R2)
		trim_R1 = sample_folder + '/' + sample_name + '_trim_R1.fastq'
		orphan_R1 = sample_folder + '/' + sample_name + '_orphan_R1.fastq'
		trim_R2 = sample_folder + '/' + sample_name + '_trim_R2.fastq'
		orphan_R2 = sample_folder + '/' + sample_name + '_orphan_R2.fastq'
	else:
		file_R1 = files[0]
		trim_R1 = sample_folder + '/' + sample_name + '_trim.fastq'

	## set command
	cmd = ""	
	if (len(files) == 2):
		cmd = "%s -jar %s PE -threads %s -trimlog %s %s %s %s %s %s %s ILLUMINACLIP:%s:2:30:10 LEADING:11 TRAILING:11 SLIDINGWINDOW:4:20 MINLEN:24 2> %s" %(java_path, trimmomatic_jar, threads, log_file, file_R1, file_R2, trim_R1, orphan_R1, trim_R2, orphan_R2, trimmomatic_adapters, trimmo_log)
	else:
		cmd = "%s -jar %s SE -threads %s -trimlog %s %s %s ILLUMINACLIP:%s:2:30:10 LEADING:11 TRAILING:11 SLIDINGWINDOW:4:20 MINLEN:24 2> %s" %(java_path, trimmomatic_jar, threads, log_file, file_R1, trim_R1, trimmomatic_adapters, trimmo_log)
	
	## system call & return
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
		

