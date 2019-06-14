 #!/usr/bin/env python3
'''
This code generates a protein annotation using PROKKA for the genomes assembled generated
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
import pandas as pd

## import my modules
from BacterialTyper import functions
from BacterialTyper import config

### print info prokka
def print_list_prokka():
	prokka_bin = config.get_exe('prokka')
	cmd = prokka_bin + " --listdb"
	functions.system_call(cmd)

######
def module_call(sequence_fasta, kingdom, path, name, threads):
	prokka_bin = config.get_exe('prokka')
	dirname = prokka_call(prokka_bin, sequence_fasta, kingdom, path, name, threads)
	return(dirname)	

######
def prokka_call(prokka_bin, sequence_fasta, kingdom, path, name, threads):
	## set parameters and options for prokka
	print ("\n+ Starting annotation for: %s\n" %name)
	outdir_name = path + '/' + name
	log_file = path + '/' + name + '.log'
	prokka = "%s --outdir %s --prefix %s --locustag %s --addgenes --addmrna --kingdom %s --cdsrnaolap --cpus %s %s 2> %s" %(prokka_bin, outdir_name, name, name, kingdom, threads, sequence_fasta, log_file)
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
