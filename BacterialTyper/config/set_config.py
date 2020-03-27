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
import pkg_resources

## import my modules
from BacterialTyper.scripts import functions
from BacterialTyper.config import extern_progs
from BacterialTyper.config import install_dependencies

################
## Software
################

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
		prog_ver = get_version(prog, p, Debug=debug)
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

##################
def decode(x):
	## this function is from ARIBA (https://github.com/sanger-pathogens/ariba)
	## give credit to them appropiately
	try:
		s = x.decode()
	except:
		return x

	return s

################
## Python
################
def print_module_comparison(module_name, message, color):
	"""
	Creates print message for a given module, version, message
	
	:param module_name: Name of the module
	:param message: Message to include in the print message: OK | FAILED | NOT FOUND
	:param color: Print message color: green | orange | red
	
	:type module_name: string
	:type message: string 
	:type color: string
	
	:returns: Print message
	"""
	print (colored("{:.<15}{:.>15}".format("Module: %s" %module_name, "[ ", message, " ]"), color))
			
#########

def get_python_packages():
	"""
	Retrieves the version of the python packages installed in the system.
	
	It retrieves the dependencies name conversion from file :file:`BacterialTyper/config/python/module_dependencies.csv`
	using function :func:`BacterialTyper.config.extern_progs.file_list` and :func:`BacterialTyper.scripts.functions.get_data`.
	For each module it retrieves the package version installed in the system using 
	:func:`BacterialTyper.config.set_config.check_package_version`.	
	
	:returns: Dictionary containing for each python module (key) the installed version (value).
	
	.. seealso:: This function relies on other ``BacterialTyper`` functions:
	
		- :func:`BacterialTyper.config.set_config.check_package_version`
		
		- :func:`BacterialTyper.config.extern_progs.file_list` 
		
		- :func:`BacterialTyper.scripts.functions.get_data`
	
	"""
	
	## get import names for packages:
	## some modules do not have the same name when install from pip and called from import
	file_module_dependecies = extern_progs.file_list("module_dependencies")
	module_dependecies = functions.get_data(file_module_dependecies, ',', 'index_col=0')

	my_packages_installed = []

	for each in my_packages_requirements:
		##	
		module_name = module_dependecies.loc[each, 'module_name_import']
		installed = check_package_version(module_name, Debug) ## check version installed in system
		my_packages_installed[each] = installed
		
	return (my_packages_installed)

##################
def check_python_packages(Debug, option_install):
	"""
	This functions checks wether the packages installed in the system fulfilled the 
	minimum version specified in the configuration folder. 
	
	It uses function :func:`BacterialTyper.config.set_config.get_python packages` to
	retrieve the version of the python packages installed in the system. Then it uses
	:func:`BacterialTyper.config.extern_progs.min_package_version` to retrieve the minimum
	version specified. It compares them using function :func:`BacterialTyper.config.set_config.check_install_module`.
	
	:param Debug:
	:param option_install:
	:type Debug: boolean
	:type option_install: string
	
	:returns: Print messages if packages are installed.
	
	.. seealso:: This function relies on other ``BacterialTyper`` functions:
	
		- :func:`BacterialTyper.config.set_config.get_python packages`
		
		- :func:`BacterialTyper.config.set_config.check_install_module`
		
		- :func:`BacterialTyper.config.extern_progs.min_package_version`
		
		- :func:`BacterialTyper.config.install_dependencies.python_package_install`
		
	"""
	
	## get python packages installed
	my_packages_installed = get_python_packages()

	## min versions for packages
	my_packages_requirements = extern_progs.min_package_version()

	## check each package
	for each in my_packages_requirements:
		## get min version	
		min_version = my_packages_requirements[each]

		## get version installed in system
		installed = my_packages_installed[each] 

		## check if installed
		message = check_install_module(installed, module_name, min_version, Debug)

		if (message == 'OK'):
			continue
		else:			
			if (option_install == 'install'):  # try to install
				installed = install_dependencies.python_package_install(module_name, min_version)
				message2 = check_install_module(installed, module_name, min_version, Debug)
				if (message2 == 'OK'):
					continue
				else:			
					print ("+ Attent to install package: ", module_name, " failed. Install it manually to continue with BacterialTyper\n\n")
			else:
				print ("+ Please install manually package: ", module_name, " to continue with BacterialTyper\n\n")

##################
def check_install_module(installed, module_name, min_version, Debug):
	"""
	
	"""
	 ## Not installed
	if (installed == 'n.a.'):
		message = 'NOT FOUND'
		color = 'red'
		print_module_comparison(module_name, message, color)
		
	# check version
	elif LooseVersion(installed) >= LooseVersion(min_version):
		message = 'OK'
		color = 'green'
		print_module_comparison(module_name, message, color)
	
	else:
		message = 'FAILED'
		color = 'orange'
		print_module_comparison(module_name, message, color)

	## return message
	return (message)
	
#########
def check_package_version(package, Debug):
	## this function is from ARIBA (https://github.com/sanger-pathogens/ariba)
	## give credit to them appropiately

	try:
		version = pkg_resources.get_distribution(package).version
		if (Debug):
			print ("1st method: pkg_resources.get_distribution(package).version")
	except:

		try:
			exec('import ' + package)
			version = eval(package + '.__version__')
			if (Debug):
				print ("2nd method: exec('import ' + package); version = eval(package + '.__version__')")

		except:
			#version = pkg_resources.resource_filename(package, 'version.py')
			#print(version)

			try:
				if (Debug):
					print ("3rd method: pkg_resources.resource_filename(package, 'version.py')")

				version = pkg_resources.resource_filename(package, 'version.py')
			except:
				version = 'n.a.'

	if (Debug):
		print ('Package:', package)
		print ('Version:', version)

	return(version)

################
## Perl
################
def check_perl_packages(Debug, option_install):
	return()

