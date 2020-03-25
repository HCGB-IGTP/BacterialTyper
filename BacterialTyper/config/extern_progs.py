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
	return(functions.get_data(dependencies_file, ',', 'index_col=0'))

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
	return(dependencies_df.loc[soft,"soft_name"])

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
	return(dependencies_df.loc[soft,"min_version"])

##################

##################
### Python packages
##################
def min_package_version():
	"""Returns a dictionary containing minimum version for each python package.

	Reads information from :file:`BacterialTyper/config/python_requirements.csv`.

	:returns: dictionary
	"""
	## read from file: prog2default.csv
	python_modules = file_list("python/python_requirements")
	package_min_versions = functions.file2dictionary(python_modules, "==")

	return(package_min_versions)
##################

def python_packages_dependencies():
	## ToDo set automatic from pip list
	python_packages_BacterialTyper = ('ariba', 'bs4', 'dendropy', 'pyfastaq', 'pymummer', 'pysam')

##################
def return_min_version_soft_package(package):
	version_package = min_package_version()
	return (version_package[package])

##################
def print_package_version():
	my_packages = min_package_version()
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


def main():

	data = read_dependencies()
	print(data)

	print_package_version()

	print (return_min_version_soft_package('ariba'))
	print (return_min_version_soft('kma'))


######
'''******************************************'''
if __name__== "__main__":
	main()

