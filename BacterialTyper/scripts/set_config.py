#!/usr/bin/env python3
##############################################################
## Jose F. Sanchez			      							##
## Copyright (C) 2019-2020 Lauro Sumoy Lab, IGTP, Spain		##
##############################################################
"""Provides configuration for the pipeline."""
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
from BacterialTyper.scripts import functions
from BacterialTyper.scripts import extern_progs
from BacterialTyper.scripts import install_dependencies

##################
def prog_to_default():
	"""Returns a dictionary containing file name for each software."""
	
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
		'mash':'mash',
		'xtract':'xtract'
	}
	return(program_to_default)

##################
def min_version_programs():
	"""Returns a dictionary containing minimum version for each software.
	
	Reads information from :file:`config/main/software_requirements.csv`.
	
	:returns: dictionary
	"""
	
	## ToDO: read requirements file:: software_requirements.csv
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
		'mash':'2.1.1',
		'efetch':'11.7',
		'esearch':'11.7',
		'xtract':'11.7',
				
		##
		'python':'3.6'
	}
	
	return min_versions

##################
def min_package_version():
	"""Returns a dictionary containing minimum version for each python package.

	Reads information from :file:`config/main/python_requirements.csv`.
	
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
def get_exe(prog, Debug=False):
	"""Return absolute path of the executable program requested.
	
	Given a program name it returns its exectuable to be called. It has to fulfilled a minimun version specified.
	
	:param prog: Software name
	
	:type prog: string
	
	:returns: Absolute path for the executable requested
	
	:warning: if no executable available in $PATH or not fulfilling the expected version.
	
	.. seealso:: This function depends on other BacterialTyper functions:
		
		- :func:`BacterialTyper.config`

		- :func:`BacterialTyper.scripts.extern_progs.return_default`
		
		- :func:`BacterialTyper.scripts.extern_progs.retrun_min_version`
		
		- :func:`BacterialTyper.scripts.extern_progs.get_version`
		
		- :func:`BacterialTyper.scripts.extern_progs.my_which`
		
	.. attention:: Be aware of Copyright
	
		The code implemented here was retrieved and modified from ARIBA (https://github.com/sanger-pathogens/ariba)
		
		Give them credit accordingly.
		
	"""
	exe = ""
	if prog in os.environ: 
		exe = os.environ[env_var] ## python environent variables
	else:
		exe = extern_progs.return_default(prog) ## install in the system

	## get paths
	exe_path_tmp = my_which(exe)
	#print (exe_path_tmp)

	## get min_version
	min_version = extern_progs.return_min_version(prog)
	#print ("Min version: ", min_version)
	
	## debugging messages
	debug=False
	if Debug:
		debug=True
	
	for p in exe_path_tmp:
		prog_ver = extern_progs.get_version(prog, p, Debug=debug)
		#print ("Path: ", p , "\nVersion: ", prog_ver)
		if (prog_ver == 'n.a.'):
			continue
		
		if LooseVersion(prog_ver) >= LooseVersion(min_version):
			return (p)
	
	if (len(exe_path_tmp) == 0):
		print(colored("\n**ERROR: Programme %s could not be found." % prog,'red'))
		exit()
	else:
		print(colored("\n**ERROR: Programme %s version smaller than minimun version expected %s." %(prog,min_version),'red'))
		exit()
			
	return('ERROR')

#################
def _access_check(fn, mode=os.F_OK | os.X_OK):
	"""Check exec permission
	
	This function checks wether a given path is a folder or file and wether it is 
	executable and accessible. It also works if a java jar file provided.
	
	:param fn: Absolute path file
	:param mode: Value to pass as the mode parameter of access()
	
	:type fn: string
	:type mode: string
	
	`mode` defaults to:
	
		- os.F_OK: Value to pass as the mode parameter of access() to test the existence of path.
		
		- os.X_OK: Value to include in the mode parameter of access() to determine if path can be executed.
	
	.. attention:: Be aware of Copyright
	
		The code implemented here was retrieved and modified from shutil (https://github.com/python/cpython/blob/master/Lib/shutil.py).
		
		Give them credit accordingly.
		
		We modified the code to work if java jar files provided.
	"""
	## the original code belongs to shutil, slightly modified here
	# https://github.com/python/cpython/blob/master/Lib/shutil.py
	
	if os.path.isdir(fn):
		return False

	if os.path.exists(fn):
		if fn.endswith('.jar'):
			return True

		if os.access(fn, mode):
			return True
	
#################
def my_which(cmd):
	"""Return the absolute path to the executable
	
	Given a command return the absolute path(s), if any.
	
	:param cmd: Software command name 
	:returns: List of absolute paths(s) of the given command.
	
	.. attention:: Be aware of Copyright
	
		The code implemented here was retrieved and modified from shutil (https://github.com/python/cpython/blob/master/Lib/shutil.py).
		
		Give them credit accordingly.
		
		We modified the code to return multiple paths in a list if available different installed binaries in $PATH.

	"""
	# If we're given a path with a directory part, look it up directly rather
	# than referring to PATH directories. This includes checking relative to the
	# current directory, e.g. ./script
	if os.path.dirname(cmd):
		if _access_check(cmd):
			return cmd
		return None

	use_bytes = isinstance(cmd, bytes)

	path=None

	if path is None:
		path = os.environ.get("PATH", None)

	if path is None:
		try:
			path = os.confstr("CS_PATH")
		except (AttributeError, ValueError):
			# os.confstr() or CS_PATH is not available
			path = os.defpath
		# bpo-35755: Don't use os.defpath if the PATH environment variable is
		# set to an empty string

	# PATH='' doesn't match, whereas PATH=':' looks in the current directory
	if not path:
		return None

	if use_bytes:
		path = os.fsencode(path)
		path = path.split(os.fsencode(os.pathsep))
	else:
		path = os.fsdecode(path)
		path = path.split(os.pathsep)

	# On other platforms you don't have things like PATHEXT to tell you
	# what file suffixes are executable, so just pass on cmd as-is.
	files = [cmd]

	return_paths = [] ## modification
	seen = set()
	for dir in path:
		normdir = os.path.normcase(dir)
		#print ("Normdir: ", normdir)
		if not normdir in seen:
			seen.add(normdir)
			for thefile in files:
				name = os.path.join(dir, thefile)
				#print ("Name: ", name)
				if _access_check(name):
					## return (name) ## previously, it would only return the first item
					return_paths.append(name) ## modification
	
	if (len(return_paths) >= 1):
		return return_paths
	else:
		return None

