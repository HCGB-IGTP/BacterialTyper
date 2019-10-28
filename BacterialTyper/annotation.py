#!/usr/bin/env python3
##########################################################
## Jose F. Sanchez										##
## Copyright (C) 2019 Lauro Sumoy Lab, IGTP, Spain		##
##########################################################
'''
Generates a protein annotation using PROKKA for the genomes assembled generated
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
from termcolor import colored

## import my modules
from BacterialTyper import functions
from BacterialTyper import config

### print info prokka
def print_list_prokka():
	prokka_bin = config.get_exe('prokka')
	cmd = prokka_bin + " --listdb"
	functions.system_call(cmd)

######
def module_call(sequence_fasta, kingdom, genus, path, name, threads):
	## check if previously assembled and succeeded
	filename_stamp = path + '/.success'

	if os.path.isdir(path):
		if os.path.isfile(filename_stamp):
			stamp =	functions.read_time_stamp(filename_stamp)
			print (colored("\tA previous command generated results on: %s" %stamp, 'yellow'))
			return ()
	
	## call prokka
	prokka_bin = config.get_exe('prokka')
	dirname = prokka_call(prokka_bin, sequence_fasta, kingdom, genus, path, name, threads)

	## success stamps
	filename_stamp = path + '/.success'
	stamp =	functions.print_time_stamp(filename_stamp)

	return(dirname)	

######
def prokka_call(prokka_bin, sequence_fasta, kingdom, genus, outdir_name, name, threads):
	## set parameters and options for prokka
	print ("\n+ Starting annotation for: %s\n" %name)
	log_file = outdir_name + '/run.log'
	options = "--cdsrnaolap --addgenes --addmrna --kingdom " + kingdom
	if genus != "Other":
		options = options + " --usegenus --genus " + genus
	prokka = "%s --force --outdir %s --prefix %s --locustag %s %s --cpus %s %s 2> %s" %(prokka_bin, outdir_name, name, name, options, threads, sequence_fasta, log_file)
	functions.system_call(prokka)
	return(outdir_name)

######
def help_options():
	print ("\nUSAGE: python %s sequence_file path name CPUs kingdom prokka_bin\n"  %os.path.realpath(__file__))

######
def main():

 	## control if options provided or help
	if len(sys.argv) > 1:
		print ("")
	else:
		help_options()
		exit()
	
	## argv
	seq_file = os.path.abspath(argv[1])
	folder = os.path.abspath(argv[2])
	name = argv[3]
	threads = int(argv[4])
	kingdom = argv[5]
	prokka_bin = argv[6]
	
	## call
	dir_annot=prokka_call(prokka_bin, seq_file, kingdom, folder, name, threads)
	print ("\n+ Annotation for sample %s has been generated in folder: %s" %(name, dir_annot))

######
'''******************************************'''
if __name__== "__main__":
	main()
