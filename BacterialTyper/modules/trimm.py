#!/usr/bin/env python3
'''
This code calls Trimmomatic for the trimming of sequence adapter within fastq reads.
Jose F. Sanchez
Copyright (C) 2019 Lauro Sumoy Lab, IGTP, Spain
'''

## import my modules
from BacterialTyper import trimmomatic_call
from BacterialTyper import multiQC_report


def run(options):
	
	
	## check folder provided and get files
	print ("Hello world!\n")
	print (options)
	print (options.input_folder)

	if (options.pair):
		print ("paired-end mode")

	else:
		print ("sinle end mode")
		
	## a folder provided containing files to trim
	
	#trimmomatic_call.trimmo_module(file_R1, file_R2, path_name, sample_name, threads)
	
	#multiQC_report.multiQC_module_call(givenList, name, path)
