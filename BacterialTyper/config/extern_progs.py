#!/usr/bin/env python3
##########################################################
## this modules is an idea from ARIBA (https://github.com/sanger-pathogens/ariba)
## give credit to them appropiately
##
## Jose F. Sanchez										##
## Copyright (C) 2019 Lauro Sumoy Lab, IGTP, Spain		##
##########################################################
"""
Provides external programs details and configuration.
"""

## useful imports
import os
import io
import sys
import re
import shutil
from io import open
from sys import argv
import subprocess
import pandas as pd
from termcolor import colored
from distutils.version import LooseVersion
import pkg_resources

## import my modules
from BacterialTyper.scripts import functions
from BacterialTyper.config import set_config
from BacterialTyper.config import install_dependencies

prog_to_version_cmd = {
	'augustus':('--version', re.compile('AUGUSTUS.*\(([0-9\.]+)\).*')),
	'ariba':('version', re.compile('ARIBA version:\s([0-9\.]+)')),
	'blastn':('-version', re.compile('blastn:\s([0-9\.]+)')),
	'busco':('--version', re.compile('BUSCO\s([0-9\.]+)')),
	'busco_plot':('--version', re.compile('BUSCO\s([0-9\.]+)')),
	'bowtie2': ('--version', re.compile('.*bowtie2.*version (.*)$')),
	'cdhit': ('', re.compile('CD-HIT version ([0-9\.]+) \(')),
	'fastqc':('-v', re.compile('FastQC\sv([0-9\.]+)')),
	'hmmsearch':('-h', re.compile('^\#.*MER\s(.*);.*')),
	'java':('-version', re.compile('version\s\"([0-9\..*\_.*]+)\"')),
	'kma':('-v', re.compile('KMA-([0-9\.]+)')),
	'prokka':('-v', re.compile('prokka\s([0-9\.]+)')),
	'makeblastdb':('-version', re.compile('makeblastdb:\s([0-9\.]+)')),
	'nucmer': ('--version', re.compile('([0-9]+\.[0-9\.]+.*$)')),
	'Rscript':('--version', re.compile('.*version\s([0-9\.]+).*')),
	'spades': ('--version', re.compile('SPAdes\s+v([0-9\.]+)')),
	'tblastn':('-version', re.compile('tblastn:\s([0-9\.]+)')),
	'multiqc':('--version', re.compile('multiqc, version\s([0-9\.]+)')),
	'trimmomatic':('-version', re.compile('([0-9\.]+)')),
	'mash':('', re.compile('Mash version ([0-9\.]+)')),
	'esearch':('-help', re.compile('esearch ([0-9\.]+)')),
	'efetch':('-help', re.compile('efetch ([0-9\.]+)')),
	'xtract':('-version', re.compile('([0-9\.]+)')),
	
}

####################################################################
def file_list(wanted_data):
	"""
	Retrieves information of additional files under folder ``BacterialTyper/config``.
	
	Using :func:`BacterialTyper.scripts.functions.get_fullpath_list` retrieves absolute
	path for file of interest.
	
	:param wanted_data: name for file
	:type wanted_data: string
	
	:returns: Absolute path for file wanted
	
	"""

	config_folder = os.path.dirname(os.path.realpath(__file__))
	listOffiles = functions.get_fullpath_list(config_folder)
	for f in listOffiles:
		name = os.path.splitext(os.path.basename(f))[0]
		if (name == wanted_data):
			return (f)

##################
## Software
##################
def read_dependencies():
	"""Returns a dictionary containing the executable name for each software.
	
	It uses :func:`BacterialTyper.config.extern_progs.file_list` to retrieve absolute path
	for file :file:`BacterialTyper/config/dependencies.csv`. It then reads csv into pandas
	dataframe using :func:`BacterialTyper.scripts.functions.get_data` and returns it.
	
	.. seealso:: This function depends on other BacterialTyper functions:
	
		- :func:`BacterialTyper.config.extern_progs.file_list`
		
		- :func:`BacterialTyper.scripts.functions.get_data`	
	"""
	
	## read from file: prog2default.csv
	dependencies_file = file_list("dependencies")
	return(functions.get_data(dependencies_file, ',', index_col=0))

	
def return_defatult_soft(soft):
	"""Returns default name for a given software name
	
	For some software we provide a shorter name for the software. Here, we read file
	:file:`BacterialTyper/config/dependencies.csv` using :func:`BacterialTyper.config.extern_progs.read_dependencies`
	and retrieve original software name.
	
	:param soft: Software name
	:type soft: string
	
	:returns: String with original software name
	
	.. seealso:: This function depends on other BacterialTyper functions:
	
		- :func:`BacterialTyper.config.extern_progs.read_dependencies`	
	
	"""
	dependencies_df = read_dependencies()
	return(dependencies_df[soft]["soft_name"])	
	
##################
def return_min_version_soft(soft):
	"""Retrieve version for a given software
	
	Retrieves minimum version for the software of interest stored in :file:`BacterialTyper/config/dependencies.csv'.
	It reads file using :func:`BacterialTyper.config.extern_progs.read_dependencies`
	and retrieve minimum version required.
	
	:param soft: Software name
	:type soft: string
	
	:returns: String with minimum version
	
	.. seealso:: This function depends on other BacterialTyper functions:
	
		- :func:`BacterialTyper.config.extern_progs.read_dependencies`	
	"""
	dependencies_df = read_dependencies()
	return(dependencies_df[soft]["min_version"])	

##################

##################
### Python packages
##################
def min_package_version():
	"""Returns a dictionary containing minimum version for each python package.

	Reads information from :file:`BacterialTyper/config/python_requirements.csv`.
	
	:returns: dictionary
	"""
	
	## ToDO: read requirements file:: python_requirements.txt
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

def python_packages_dependencies():
	## ToDo set automatic from pip list
	python_packages_BacterialTyper = ('ariba', 'bs4', 'dendropy', 'pyfastaq', 'pymummer', 'pysam')

##################
def return_min_version_soft_package(package):
	version_package = config.min_package_version()
	return (version_package[package])

##################
def print_package_version():
	my_packages = config.min_package_version()
	for each in my_packages:
		print ("{:.<15}{:.>15}".format("Module: %s" %each, my_packages[each]))






##################
def print_dependencies():
	progs = {}
	prog_to_default = config.prog_to_default()
	for prog in prog_to_default:
		#print (prog)
		prog_exe = config.get_exe(prog)
		#print (prog + '\t' + prog_exe)
		prog_ver = get_version(prog, prog_exe)
		progs[prog] = [prog_exe, prog_ver]

	df_programs = pd.DataFrame.from_dict(progs, orient='index', columns=('Executable path', 'Version'))
	df_programs = df_programs.stack().str.lstrip().unstack()
	pd.set_option('display.max_colwidth', -1)
	pd.set_option('display.max_columns', None)
	print (df_programs)

