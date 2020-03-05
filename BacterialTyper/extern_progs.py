#!/usr/bin/env python3
##########################################################
## this modules is an idea from ARIBA (https://github.com/sanger-pathogens/ariba)
## give credit to them appropiately
##
## Jose F. Sanchez										##
## Copyright (C) 2019 Lauro Sumoy Lab, IGTP, Spain		##
##########################################################
'''
Provides external programs details
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
import pandas as pd
from termcolor import colored
from distutils.version import LooseVersion
import pkg_resources

## import my modules
from BacterialTyper import functions
from BacterialTyper import config
from BacterialTyper import install_dependencies

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

##################
def return_default(soft):
	dict_programs = config.prog_to_default()
	return (dict_programs[soft])
	
##################
def return_min_version(soft):
	"""Retrieve version for a given software
	
	Retrieves minimun version for the software of interest stored in :file:`config.main.software_requirements.csv' using the function :func:`BacterialTyper.config.min_version_programs`.
	
	.. seealso:: Additional information on BacterialTyper configuration and requirements
	
		- :doc:`Configuration <../../../user_guide/installing>` 

	"""
	version_programs = config.min_version_programs()
	return (version_programs[soft])

##################
def python_packages_dependencies():
	## ToDo set automatic from pip list
	python_packages_BacterialTyper = ('ariba', 'bs4', 'dendropy', 'pyfastaq', 'pymummer', 'pysam')

##################
def return_min_version_package(package):
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

##################
def decode(x):
	## this function is from ARIBA (https://github.com/sanger-pathogens/ariba)
	## give credit to them appropiately
	try:
		s = x.decode()
	except:
		return x
	
	return s

##################
def get_version(prog, path, Debug=False):
	## this function is from ARIBA (https://github.com/sanger-pathogens/ariba)
	## give credit to them appropiately
	"""Get version of software
	
	Given a program name and expected path, tries to determine its version.
	
	:param prog:
	:param path:
	:param Debug:
	
	:type prog:
	:type path:
	:type Debug:
	
	:returns: tuple (bool, string). First element True if found version ok.
	Second element is version. Returns NA message if no found and raises attention error message.
	
	.. attention:: Be aware of Copyright
	
		The code implemented here was retrieved and modified from ARIBA (https://github.com/sanger-pathogens/ariba)
		
		Give them credit accordingly.
	"""

	assert prog in prog_to_version_cmd
	args, regex = prog_to_version_cmd[prog]
	cmd = path + ' ' + args
	if prog == 'spades':
		cmd_output = subprocess.Popen(['python3', path, args], shell=False, stdout=subprocess.PIPE,stderr=subprocess.PIPE).communicate()
	elif prog == 'trimmomatic':
		java_bin = config.get_exe('java')
		java_jar = java_bin + ' -jar ' + path + ' ' + args
		cmd_output = subprocess.Popen(java_jar, shell=True, stdout=subprocess.PIPE,stderr=subprocess.PIPE).communicate()
	else:
		cmd_output = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()

	## decode command
	cmd_output = decode(cmd_output[0]).split('\n')[:-1] + decode(cmd_output[1]).split('\n')[:-1]
	
	## retrieve version information
	for line in cmd_output:
		hits = regex.search(line)
		if hits:
			return hits.group(1)
	
	if Debug:
		print (colored('Attention: I tried to get the version of ' + prog + ' with: "' + cmd + '" and the output didn\'t match this regular expression: "' + regex.pattern + '"', 'red'))

	return("n.a.")


#########
def check_python_packages(Debug, install):
	## get all packages and min versions
	my_packages = config.min_package_version()
	for each in my_packages:
		##	
		min_version = my_packages[each]
		installed = check_package_version(each) ## check version installed in system

		## Not installed
		if (installed == 'n.a.'):
			print (colored("{:.<15}{:.>15}".format("Module: %s" %each, "[ NOT FOUND ]"), 'red'))
			if (Debug):
				print ("\n**", each, min_version, installed, " **")			
			if (install): # try to install
				installed = install_dependencies.python_package_install(each, min_version)
				if (Debug):
					print ("\n**", each, min_version, installed, " **")			
			else:
				continue

		# check version
		if LooseVersion(installed) >= LooseVersion(min_version):
			print (colored("{:.<15}{:.>15}".format("Module: %s" %each, "[ OK ]"), 'green'))

		else:
			print (colored("{:.<15}{:.>15}".format("Module: %s" %each, "[ FAILED ]"), 'red'))
			#print (colored("Package %s\t[ FAILED ]" % each,'red'))
			if (install):  # try to install
				installed = install_dependencies.python_package_install(each, min_version)
				if (Debug):
					print ("\n**", each, min_version, installed, " **")	
				if LooseVersion(installed) >= LooseVersion(min_version):
					print (colored("{:.<15}{:.>15}".format("Module: %s" %each, "[ OK ]"), 'green'))
				else:
					print (colored("{:.<15}{:.>15}".format("Module: %s" %each, "[ FAILED (II) ]"), 'red'))
					#print (colored("Package %s\t[ FAILED (II) ]" % each,'red'))
					print ("+ Please install manually package: ", each, "\n\n")
			else:
				continue

#########
def check_package_version(package):
	## this function is from ARIBA (https://github.com/sanger-pathogens/ariba)
	## give credit to them appropiately

	try:
		version = pkg_resources.get_distribution(package).version
		return (version)
	except:
		
		try:
			exec('import ' + package)
			version = eval(package + '.__version__')
			return (version)

		except:
			return ('n.a.')
		

