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

## import my modules
from BacterialTyper import functions
from BacterialTyper import config

prog_to_default = {
	'ariba':'ariba',
   	'bowtie2': 'bowtie2',
   	'cdhit': 'cd-hit-est',
   	'nucmer' : 'nucmer',
   	'spades' : 'spades.py',
   	'kma':'kma',
   	'fastqc':'fastqc'
   	
   	##	blastn
   	##	makeblastdb
	##	bowtie2
	##	BUSCO
	##	augustus
	##	prokka
	##	trimmomatic
	
	## plasmid id
	##	bedtools
	##	samtools
	##	circos
	##	plasmidID

}
	
	
prog_to_version_cmd = {
	'bowtie2': ('--version', re.compile('.*bowtie2.*version (.*)$')),
	'cdhit': ('', re.compile('CD-HIT version ([0-9\.]+) \(')),
	'nucmer': ('--version', re.compile('([0-9]+\.[0-9\.]+.*$)')),
	'spades': ('--version', re.compile('SPAdes\s+v([0-9\.]+)')),
	'ariba':('version', re.compile('ARIBA version:\s([0-9\.]+)')),
	'kma':('-v', re.compile('KMA-([0-9\.]+)')),
	'fastqc':('-v', re.compile('FastQC\sv([0-9\.]+)'))
}


min_versions = {
	'bowtie2': '2.1.0',
	'cdhit': '4.6',
	'nucmer': '3.1',
	'spades': '3.11.0',
	'kma':'1.2.2',
	'fastqc':'0.11.4'
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
	for prog in prog_to_default:
		prog_exe = get_exe(prog)
		prog_ver = get_version(prog, prog_exe)
		progs[prog] = [prog_exe, prog_ver]

	df_programs = pd.DataFrame.from_dict(progs, orient='index')
	df_programs = df_programs.stack().str.lstrip().unstack()
	print (df_programs)

##################
def get_exe(prog):
	## this function is from ARIBA (https://github.com/sanger-pathogens/ariba)
	## give credit to them appropiately
	'''Given a program name, return what we expect its exectuable to be called'''
	exe = ""
	if prog in os.environ: 
		exe = os.environ[env_var] ## python environent variables
	else:
		exe = prog_to_default[prog] ## install in the system

	return(shutil.which(exe)) ## return which path


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
	else:
		cmd_output = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()

	cmd_output = decode(cmd_output[0]).split('\n')[:-1] + decode(cmd_output[1]).split('\n')[:-1]

	for line in cmd_output:
		hits = regex.search(line)
		if hits:
			return hits.group(1)
	
	return 'ERROR - I tried to get the version of ' + prog + ' with: "' + cmd + '" and the output didn\'t match this regular expression: "' + regex.pattern + '"'
	

