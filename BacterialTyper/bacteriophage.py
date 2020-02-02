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

## import my modules
from BacterialTyper import functions
from BacterialTyper import config
import PhiSpy
from PhiSpy_tools import genbank_to_seed

perlDir = os.path.dirname(os.path.realpath(__file__)) + '/other_tools/perl'
sys.path.append(perlDir)
get_longContigs_script = perlDir + '/get-long-contigs.pl'

######
def get_long_seqs(sequence_file, output_file):
	cmd="perl %s %s %s > %s" %(get_longContigs_script, sequence_file, 2000, output_file)
	return(functions.system_call(cmd))

######
def convert_Genbank2seed():
	return()

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


