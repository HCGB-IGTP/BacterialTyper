#!/usr/bin/env python3
##########################################################
## Jose F. Sanchez									  ##
## Copyright (C) 2019 Lauro Sumoy Lab, IGTP, Spain	  ##
##########################################################
'''
Calls FASTQC analysis and parses results generated.
'''
## useful imports
import time
import io
import os
import re
import sys
from sys import argv
from io import open
import pandas as pd

import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
from pandas.plotting import table
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np

## import my modules
import HCGB
import HCGB.functions.system_call_functions as HCGB_sys
from BacterialTyper.config import set_config

############
def	help_options():
	print ("\nUSAGE:\npython %s folder file1 [file2] name fastqc threads\n"  %os.path.realpath(__file__))

############
def call_fastqc(path, file1, file2, sample, fastqc_bin, threads):	
	## call system for fastqc sample given
	#name = functions.create_subfolder(sample, path)
	logFile = path + '/' + sample + '.log'
	
	if os.path.isfile(logFile):
		return ('OK')
	
	if not file2:
		cmd_fastqc = '%s --extract -t %s -o %s %s > %s 2> %s' %(fastqc_bin, threads, path, file1, logFile, logFile)
	else:
        ##print ("+ Calling fastqc for samples...")	
	       cmd_fastqc = '%s --extract -t %s -o %s %s %s > %s 2> %s' %(fastqc_bin, threads, path, file1, file2, logFile, logFile)
	
    ## send command	
	return (HCGB_sys.system_call( cmd_fastqc ))
		
	
############
def run_module_fastqc(path, files, sample, threads):	
	## Arguments provided via ARGVs
	fastqc_bin = set_config.get_exe('fastqc')

	if (len(files) == 1):
		codeReturn = call_fastqc(path, files[0], "", sample, fastqc_bin, threads)
	elif (len(files) == 2):
		codeReturn = call_fastqc(path, files[0], files[1], sample, fastqc_bin, threads)
	else:
		print("ERROR: Some error ocurred while executing fastqc. Check samples provided")
	
	if codeReturn == 'FAIL':
		exit()
	path_to_sample = path + '/' + sample
	return path_to_sample



############
def main():
	## this code runs when call as a single script

  	## control if options provided or help
	if len(sys.argv) > 1:
		print ("")
	else:
		help_options()
		exit()		

	path = os.path.abspath(argv[1])
	file1 = os.path.abspath(argv[2])
	file2 = os.path.abspath(argv[3])
	sample = argv[4]
	fastqc_bin = argv[5]
	threads = argv[5]

	## check if paired end
	if not (file2):
		print ('+ No implementation yet for single-end. Sorry.')
		exit()

	##
	path_to_sample = call_fastqc(path, file1, file2, sample, fastqc_bin, threads)
		
############
'''******************************************'''
if __name__== "__main__":
	main()







