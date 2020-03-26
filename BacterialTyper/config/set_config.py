#!/usr/bin/env python3
##############################################################
## Jose F. Sanchez			      							##
## Copyright (C) 2019-2020 Lauro Sumoy Lab, IGTP, Spain		##
##############################################################
"""Provides configuration for the pipeline.

.. seealso:: Additional information on BacterialTyper configuration and requirements
	
	- :doc:`Configuration <../../user_guide/installation/installing>` 
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
from termcolor import colored
from distutils.version import LooseVersion

## import my modules
from BacterialTyper.scripts import functions
from BacterialTyper.config import extern_progs
from BacterialTyper.config import install_dependencies

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

		- :func:`BacterialTyper.scripts.extern_progs.return_defatult_soft`
		
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
		exe = extern_progs.return_defatult_soft(prog) ## install in the system

	## get paths
	exe_path_tmp = my_which(exe)
	#print (exe_path_tmp)

	## get min_version
	min_version = extern_progs.return_min_version_soft(prog)
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

##################
def get_version(prog, path, Debug=False):
	## this function is from ARIBA (https://github.com/sanger-pathogens/ariba)
	## give credit to them appropiately
	"""Get version of software
	
	Given a program name and expected path, tries to determine its version.
	
	:param prog: Program name
	:param path: Absolute path
	:param Debug: True/False
	
	:type prog: string
	:type path: string 
	:type Debug: bool
	
	:returns: tuple (bool, string). First element True if found version ok. Second element is version. Returns NA message if no found and raises attention error message.
	
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
def check_python_packages(Debug, option_install):
	
	## min versions for packages
	my_packages = extern_progs.min_package_version()
	
	## get import nades for packages:
	## some modules do not have the same name when install from pip and called from import
	file_module_dependecies = extern_progs.file_list("module_dependencies")
	module_dependecies = functions.get_data(file_module_dependecies, ',', 'index_col=0')
		
	for each in my_packages:
		##	
		min_version = my_packages[each]
		module_name = module_dependecies.loc[each, 'module_name_import']
		installed = check_package_version(module_name) ## check version installed in system

		## Not installed
		if (installed == 'n.a.'):
			print (colored("{:.<15}{:.>15}".format("Module: %s" %module_name, "[ NOT FOUND ]"), 'red'))
			if (Debug):
				print ("\n**", module_name, min_version, installed, " **")			
			if (option_install == "install"): # try to install
				installed = install_dependencies.python_package_install(module_name, min_version)
				if (Debug):
					print ("\n**", module_name, min_version, installed, " **")			
			else:
				continue

		# check version
		if LooseVersion(installed) >= LooseVersion(min_version):
			print (colored("{:.<15}{:.>15}".format("Module: %s" %module_name, "[ OK ]"), 'green'))

		else:
			print (colored("{:.<15}{:.>15}".format("Module: %s" %module_name, "[ FAILED ]"), 'red'))
			#print (colored("Package %s\t[ FAILED ]" % module_name,'red'))
			if (option_install == 'install'):  # try to install
				installed = install_dependencies.python_package_install(module_name, min_version)
				if (Debug):
					print ("\n**", module_name, min_version, installed, " **")	
				if LooseVersion(installed) >= LooseVersion(min_version):
					print (colored("{:.<15}{:.>15}".format("Module: %s" %module_name, "[ OK ]"), 'green'))
				else:
					print (colored("{:.<15}{:.>15}".format("Module: %s" %module_name, "[ FAILED (II) ]"), 'red'))
					#print (colored("Package %s\t[ FAILED (II) ]" % module_name,'red'))
					print ("+ Please install manually package: ", module_name, "\n\n")
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
		
##################
def decode(x):
	## this function is from ARIBA (https://github.com/sanger-pathogens/ariba)
	## give credit to them appropiately
	try:
		s = x.decode()
	except:
		return x
	
	return s



