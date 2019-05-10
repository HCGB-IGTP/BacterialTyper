#usr/bin/env python
'''
This module provides configuration for the pipeline
Jose F. Sanchez
Copyright (C) 2019 Lauro Sumoy Lab, IGTP, Spain
'''
## useful imports
import os
import io
import sys
import re
import shutil
from io import open
from sys import argv
import subprocess
from termcolor import colored

## import my modules
from BacterialTyper import functions
from BacterialTyper import extern_progs
from BacterialTyper import install_dependencies

##################
def prog_to_default():

	program_to_default = {
		'ariba':'ariba',
		'bowtie2': 'bowtie2',
		'cdhit': 'cd-hit-est',
		'nucmer' : 'nucmer',
		'spades' : 'spades.py',
		'kma':'kma',
		'fastqc':'fastqc',
		'busco':'run_BUSCO.py',
		'tblastn':'tblastn',
		'blastn':'blastn',
		'makeblastdb':'makeblastdb',
		'bowtie2':'bowtie2',
		'busco':'run_BUSCO.py',
		'prokka':'prokka',

		#'trimmomatic':'trimmomatic.jar',
		'hmmsearch':'hmmsearch',

		'augustus':'augustus',
		'Rscript':'Rscript',
		'java':'java'
		
		## plasmid id
		##	bedtools
		##	samtools
		##	circos
		##	plasmidID

	}
	return(program_to_default)


##################
def return_default(soft):
	dict_programs = prog_to_default()
	return (dict_programs[soft])	

##################
def get_exe(prog):
	## this function is from ARIBA (https://github.com/sanger-pathogens/ariba)
	## give credit to them appropiately
	'''Given a program name, return what we expect its exectuable to be called'''
	exe = ""
	if prog in os.environ: 
		exe = os.environ[env_var] ## python environent variables
	else:
		exe = return_default(prog) ## install in the system

	exe_path = shutil.which(exe)
	if (exe_path):
		return(exe_path) ## return which path
	else:
		print(colored("**ERROR: Programme %s could not be found." % prog,'red'))
		return('ERROR')

