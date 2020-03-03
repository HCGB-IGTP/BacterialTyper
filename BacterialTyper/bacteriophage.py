#!/usr/bin/env python3
'''
This code generates a bacteriophage identification profile for each sample
Jose F. Sanchez
Copyright (C) 2019 Lauro Sumoy Lab, IGTP, Spain
'''
## useful imports
import time
import io
import os
import re
import sys
from sys import argv
from io import open
from termcolor import colored

## import my modules
from BacterialTyper import functions
from BacterialTyper import config

## import phispy modules
import PhiSpyModules
#from PhiSpy_tools import genbank_to_seed
# Contribution: https://github.com/linsalrob/PhiSpy/PhiSpyModules/writers.py

## [TODO]

######
def ident_bacteriophage(gbk_file, outdir, training_set, window_size, phage_genes, Debug):
	"""Call PhiSpy to identify putative bacteriophages.
	"""

	if os.path.isfile(filename_stamp):
		stamp =	functions.read_time_stamp(filename_stamp)
		print (colored("\tA previous command generated results on: %s" %stamp, 'yellow'))
	else:
	
		## debug message
		if (Debug):
			print (colored("**DEBUG: convert_Genbank2seed**", 'yellow'))
			print (" convert_Genbank2seed(gbk_file, outdir)")
			print (" convert_Genbank2seed " + "\t" + gbk_file + "\t" + outdir_gbk)

		# Call genbank conversion
		code_returned = convert_Genbank2seed(gbk_file, outdir_gbk)

	if (code_returned == 'OK'):
		#call phispy
		print ("+ Identify putative phage regions...")
	
	else:
		return('FAIL')


	## prophage_coordiantes.csv_
	## "#prophage_ID", "Contig","Start","End","attL_Start","attL_End","attR_Start","attR_End","attL_Seq","attR_Seq","Longest_Repeat_flanking_phage" 

	## prophage.gff3
	# Contribution: https://github.com/linsalrob/PhiSpy/PhiSpyModules/writers.py





######
def help_PhiSpy():
	print ("\n** PhiSpy additional information **")
	print ("phiSpy is a program for identifying prophages from among microbial genome sequences\n")
	print ("(c) 2008-2018 Sajia Akhter, Katelyn McNair, Rob Edwards, San Diego State University, San Diego, CA")
	print ("https://github.com/linsalrob/PhiSpy\n")
	
	print ("Choose among different training sets for a better phage identification:\n")
	PhiSpyModules.print_list()	
	print ("\n")

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
	
	PhiSpyModules.print_list()	
	genbank_to_seed.usage()

	
	## argv
	#seq_file = os.path.abspath(argv[1])
	#folder = os.path.abspath(argv[2])
	#name = argv[3]
	#threads = int(argv[4])
	
	## call
	#print ("")
	#folder_name = functions.create_subfolder(folder, name)
	#get_long_seqs(seq_file, folder_name + '/tmp.fasta')	
	
	#print ("\n+ Annotation of putative phages for sample %s has been generated in folder: %s" %(name, dir_annot))

######
'''******************************************'''
if __name__== "__main__":
	main()
