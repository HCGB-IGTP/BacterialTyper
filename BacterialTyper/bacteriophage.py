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
import PhiSpy
from PhiSpy_tools import genbank_to_seed

perlDir = os.path.dirname(os.path.realpath(__file__)) + '/other_tools/perl'
sys.path.append(perlDir)
get_longContigs_script = perlDir + '/get-long-contigs.pl'

## [TODO]

######
def get_long_seqs(sequence_file, output_file):
	cmd="perl %s %s %s > %s" %(get_longContigs_script, sequence_file, 2000, output_file)
	return(functions.system_call(cmd))

######
def convert_Genbank2seed(gbk_file, path):
	print (colored("\n\n***** TODO: Debug Implement genbank2seed call *****\n\n", 'red'))	
	code = genbank_to_seed.convert_contigs(["", gbk_file, path])
	if (code == 'OK'):
		## success stamps
		filename_stamp = path + '/.success_genbank2seed'
		stamp =	functions.print_time_stamp(filename_stamp)

	return(code)

######
def ident_bacteriophage(gbk_file, outdir, training_set, window_size, phage_genes, Debug):
	
	## convert genbank2seed
	print ("+ Convert genbank file into seed format")
	
	outdir_gbk = outdir + '/gbk2seed'
	
	## check if previously assembled and succeeded
	filename_stamp = outdir_gbk + '/.success_genbank2seed'

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

######
def help_PhiSpy():
	print ("\n** phiSpy additional information **")
	print ("phiSpy is a program for identifying prophages from among microbial genome sequences\n")
	print ("(c) 2008-2018 Sajia Akhter, Katelyn McNair, Rob Edwards, San Diego State University, San Diego, CA")
	print ("https://github.com/linsalrob/PhiSpy\n")
	
	print ("Choose among different training sets for a better phage identification:\n")
	PhiSpy.print_list()	
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
	
	PhiSpy.print_list()	
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


