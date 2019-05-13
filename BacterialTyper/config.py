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
from distutils.version import LooseVersion

## import my modules
from BacterialTyper import functions
from BacterialTyper import extern_progs
from BacterialTyper import install_dependencies

##################
def prog_to_default():

	program_to_default = {
		'ariba':'ariba',
		'augustus':'augustus',
		'blastn':'blastn',
		'bowtie2': 'bowtie2',
		'busco':'run_BUSCO.py',
		'cdhit': 'cd-hit-est',
		'fastqc':'fastqc',
		'hmmsearch':'hmmsearch',
		'java':'java',
		'kma':'kma',
		'prokka':'prokka',
		'makeblastdb':'makeblastdb',
		'nucmer' : 'nucmer',
		'Rscript':'Rscript',
		'spades' : 'spades.py',
		'tblastn':'tblastn',
		'trimmomatic':'trimmomatic.jar'
	}
	return(program_to_default)

	## plasmid id
	##	bedtools
	##	samtools
	##	circos
	##	plasmidID


##################
def return_default(soft):
	dict_programs = prog_to_default()
	return (dict_programs[soft])

##################
def min_version_programs():

	min_versions = { ## update
		'ariba':'2.13.5',
		'augustus':'3.2.1',		
		'blastn':'2.5',
		'bowtie2': '2.1.0',
		'busco':'3.1.0',
		'cdhit': '4.6',
		'fastqc':'0.11.4',
		'hmmsearch':'3.1b2',
		'java':'1.8.0_172',
		'kma':'1.2.2',
		'prokka':'1.12',
		'makeblastdb':'2.5',
		'nucmer': '3.1',
		'Rscript':'3.5.1',
		'spades':'3.9.0',		
		'tblastn':'2.5',
		'trimmomatic':'0.36'
	}
	
	return min_versions


##################
def return_min_version(soft):
	version_programs = min_version_programs()
	return (version_programs[soft])

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

	## get paths
	exe_path_tmp = functions.my_which(exe)
	#print (exe_path_tmp)

	## get min_version
	min_version = return_min_version(prog)
	#print ("Min version: ", min_version)
	
	for p in exe_path_tmp:
		prog_ver = extern_progs.get_version(prog, p)
		#print ("Path: ", p , "\nVersion: ", prog_ver)
		if LooseVersion(prog_ver) >= LooseVersion(min_version):
			return (p)
	
	if (len(exe_path_tmp) == 0):
		print(colored("**ERROR: Programme %s could not be found." % prog,'red'))
	else:
		print(colored("**ERROR: Programme %s version smaller than minimun version expected %s." %(prog,min_version),'red'))
	
	return('ERROR')

