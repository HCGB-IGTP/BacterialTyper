#usr/bin/env python
'''
This module provides external programs details
Jose F. Sanchez
Copyright (C) 2019 Lauro Sumoy Lab, IGTP, Spain
'''
## this modules is an idea from ARIBA (https://github.com/sanger-pathogens/ariba)
## give credit to them appropiately

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

## import my modules
from BacterialTyper import functions
from BacterialTyper import config

prog_to_version_cmd = {
	'augustus':('--version', re.compile('AUGUSTUS.*\(([0-9\.]+)\).*')),
	'ariba':('version', re.compile('ARIBA version:\s([0-9\.]+)')),
	'blastn':('-version', re.compile('blastn:\s([0-9\.]+)')),
	'busco':('--version', re.compile('BUSCO\s([0-9\.]+)')),
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
	'trimmomatic':('-version', re.compile('([0-9\.]+)'))
}

package_min_versions = {
    'bs4': '4.1.0',
    'dendropy': '4.1.0',
    'pyfastaq': '3.12.0',
    'pysam': '0.8.1',
    'pymummer' : '0.7.1',
}


##################
def dependencies():
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
def get_version(prog, path):
	## this function is from ARIBA (https://github.com/sanger-pathogens/ariba)
	## give credit to them appropiately
	'''Given a program name and expected path, tries to determine its version.
	Returns tuple (bool, version). First element True iff found version ok.
	Second element is version string (if found), otherwise an error message'''
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
	
	print (colored('ERROR - I tried to get the version of ' + prog + ' with: "' + cmd + '" and the output didn\'t match this regular expression: "' + regex.pattern + '"', 'red'))
	return("n.a.")


	

