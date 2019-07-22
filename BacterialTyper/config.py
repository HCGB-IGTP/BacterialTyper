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
		'busco_plot':'generate_plot.py',
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
		'multiqc':'multiqc',
		'trimmomatic':'trimmomatic.jar',
		'efetch':'efetch',
		'esearch':'esearch',
		'xtract':'xtract'
	}
	return(program_to_default)

##################
def min_version_programs():

	min_versions = { ## update
		'ariba':'2.13.5',
		'augustus':'3.2.1',		
		'blastn':'2.5',
		'bowtie2': '2.1.0',
		'busco':'3.1.0',
		'busco_plot':'3.1.0',
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
		'trimmomatic':'0.36',
		'multiqc':'1.7',
		
		'efetch':'11.7',
		'esearch':'11.7',
		'xtract':'11.7',
				
		##
		'python':'3.5'
	}
	
	return min_versions

##################
def min_package_version():
	package_min_versions = {
		'appdirs':'1.4.3',
		'ariba':'2.13.5',
		'Bio':'1.73', ## biopython
		'bs4':'4.7.1', #beautifulsoup4
		'certifi':'2019.3.9',
		'chardet':'3.0.4',
		'click':'7.0',
		'colormath':'3.0.0',
		'configparser':'3.7.4',
		'cycler':'0.10.0',
		'cython':'0.29.6',
		'decorator':'4.4.0',
		'dendropy':'4.4.0',
		'et_xmlfile':'1.0.1',
		'ete3':'3.1.1',
		'fastqcparser':'1.1',
		'filehash':'0.1.dev3',
		'future':'0.17.1',
		'idna':'2.8',
		'jdcal':'1.4.1',
		'jinja2':'2.10.1',
		'kiwisolver':'1.0.1',
		'lzstring':'1.0.4',
		'markdown':'3.1',
		'markupsafe':'1.1.1',
		'matplotlib':'2.2.4',
		'multiqc':'1.7',
		'ncbi_genome_download':'0.2.9',
		'networkx':'2.2',
		'numpy':'1.16.2',
		'openpyxl':'2.6.2',
		'pandas':'0.24.2',
		'patoolib':'1.12',
		'pyfastaq':'3.17.0',
		'pymummer':'0.10.3',
		'pyparsing':'2.4.0',
		'pysam':'0.15.2',
		'dateutil':'2.8.0',
		'python-magic':'0.4.15',
		'pytz':'2018.9',
		'yaml':'5.1', #pyyaml
		'requests':'2.21.0',
		'scipy':'1.2.1',
		'simplejson':'3.16.0',
		'six':'1.12.0',
		'soupsieve':'1.9',
		'spectra':'0.0.11',
		'termcolor':'1.1.0',
		'urllib3':'1.24.1',
		'wget':'3.2',
		'xlrd':'1.2.0',
		'xlsxwriter':'1.1.7',
		'xlwt':'1.3.0'
	}
	return package_min_versions

##################
def get_exe(prog):
	## this function is from ARIBA (https://github.com/sanger-pathogens/ariba)
	## give credit to them appropiately
	'''Given a program name, return what we expect its exectuable to be called'''
	exe = ""
	if prog in os.environ: 
		exe = os.environ[env_var] ## python environent variables
	else:
		exe = extern_progs.return_default(prog) ## install in the system

	## get paths
	exe_path_tmp = functions.my_which(exe)
	#print (exe_path_tmp)

	## get min_version
	min_version = extern_progs.return_min_version(prog)
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

