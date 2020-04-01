#!/usr/bin/env python3
##########################################################
## Jose F. Sanchez					##
## Copyright (C) 2019-2020 Lauro Sumoy Lab, IGTP, Spain	##
##########################################################
"""
Provides configuration for the pipeline.

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
def get_exe(prog, Debug=False, Return_Version=False):
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

	## debug message
	if (Debug):
		print(colored("** Debug: exe: %s" %exe,'red'))
		print(colored("** Debug: exe_path_tmp: %s" %exe_path_tmp,'red'))

	## get min_version
	min_version = extern_progs.return_min_version_soft(prog)
	## debug message
	if (Debug):
		print(colored("** Debug: min_version: %s" %min_version,'red'))

	## return if not available
	if not exe_path_tmp:
		if (Return_Version):
			return ('ERROR', 'n.a.')
		else:
			return('ERROR')

	## Loop for all possibilities
	for p in exe_path_tmp:
		prog_ver = get_version(prog, p, Debug=Debug)

		if (Debug):
			print (colored("** Debug: Sotf: %s\nPath: %s\nVersion: %s" %(prog, p, prog_ver), 'red'))

		if (prog_ver == 'n.a.'):
			continue

		if LooseVersion(prog_ver) >= LooseVersion(min_version):
			if (Return_Version):
				return (p, prog_ver)
			else:
				return (p)

	if (len(exe_path_tmp) == 0):
		print(colored("\n**ERROR: Programme %s could not be found." % prog,'red'))
		exit()
	else:
		print(colored("\n**ERROR: Programme %s version smaller than minimun version expected %s." %(prog,min_version),'red'))
		exit()

	if (Return_Version):
		return ('ERROR', 'n.a.')
	else:
		return('ERROR')

################
def access_check(fn, mode=os.F_OK | os.X_OK):
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

	#if os.path.isdir(fn):
	#	return False

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
		if access_check(cmd):
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
				if access_check(name):
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

	## read dependencies information
	dependencies_pd = extern_progs.read_dependencies()

	## get information for prog
	regex = re.compile(dependencies_pd.loc[prog, 'get_version'])
	args = dependencies_pd.loc[prog, 'version_cmd']
	cmd = path + ' ' + args

	## debug messages
	if (Debug):
		print(colored("** Debug: regex: %s" %regex,'red'))
		print(colored("** Debug: args: %s" %args, 'red'))

	if prog == 'spades':
		cmd_output = subprocess.Popen(['python3', path, args], shell=False, stdout=subprocess.PIPE,stderr=subprocess.PIPE).communicate()
	elif prog == 'trimmomatic':
		java_bin = set_config.get_exe('java')
		java_jar = java_bin + ' -jar ' + path + ' ' + args
		cmd_output = subprocess.Popen(java_jar, shell=True, stdout=subprocess.PIPE,stderr=subprocess.PIPE).communicate()
	else:
		cmd_output = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()

	## decode command
	cmd_output = functions.decode(cmd_output[0]).split('\n')[:-1] + functions.decode(cmd_output[1]).split('\n')[:-1]

	## debug messages
	if (Debug):
		print(colored("** Debug: cmd_output:\n %s" %cmd_output,'red'))

	## retrieve version information
	for line in cmd_output:
		hits = regex.search(line)
		if (Debug):
			print (hits)
		if hits:
			return hits.group(1)

	if Debug:
		print (colored('Attention: I tried to get the version of ' + prog + ' with: "' + cmd + '" and the output didn\'t match this regular expression: "' + regex.pattern + '"', 'red'))

	return("n.a.")

def check_dependencies(install_option, install_path, Debug):
	"""
	"""
	
	## read dependencies information
	dependencies_pd = extern_progs.read_dependencies()

	for soft, row in dependencies_pd.iterrows():
		(soft_path, installed) = get_exe(soft, Debug=Debug, Return_Version=True)
		soft_name = row['soft_name']
		min_version = row['min_version']
		
		## debug messages
		if (Debug):
			print ("Software:", soft)
			print ("Soft name:", soft_name)
			print ("Min_Version:", min_version)
			
			print ("Soft Path: ", soft_path)
			print ("Version installed:", installed)
			
		## check if installed
		message = check_install_module(installed, soft_name, min_version, 'Software')

		if (message == 'OK'):
			continue
		else:
			if (install_option):
				if (Debug):
						print ("Install module: ", each)
	
				installed = install_dependencies.install(soft, min_version)
				message2 = check_install_module(installed, soft_name, min_version, 'Software')
				if (message2 == 'OK'):
					continue
				else:
					print ("+ Attent to install software: ", soft_name, " failed. Install it manually to continue with BacterialTyper\n\n")
			else:
				print ("+ Please install manually software: ", soft_name, " to continue with BacterialTyper\n\n")
	
################
## Python
################
def get_python_packages(Debug):
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

		- :func:`BacterialTyper.scripts.functions.file2dictionary`

	"""

	## get import names for packages:
	## some modules do not have the same name when install from pip and called from import
	file_module_dependecies = extern_progs.file_list("module_dependencies")
	module_dependencies = functions.file2dictionary(file_module_dependecies, ',')

	my_packages_installed = {}

	for each in module_dependencies:
		module_name = module_dependencies[each]
		installed = check_package_version(module_name, Debug) ## check version installed in system
		my_packages_installed[each] = installed

	return (my_packages_installed)

##################
def check_python_packages(Debug, option_install, install_path):
	"""
	This functions checks wether the packages installed in the system fulfilled the 
	minimum version specified in the configuration folder. 

	It uses function :func:`BacterialTyper.config.set_config.get_python packages` to
	retrieve the version of the python packages installed in the system. Then it uses
	:func:`BacterialTyper.config.extern_progs.min_python_module_version` to retrieve the minimum
	version specified. It compares them using function :func:`BacterialTyper.config.set_config.check_install_module`.

	:param Debug: True/False for debugging messages
	:param option_install: True/False for installing missing dependencies
	:type Debug: boolean
	:type option_install: string

	:returns: Print messages if packages are installed.

	.. seealso:: This function relies on other ``BacterialTyper`` functions:

		- :func:`BacterialTyper.config.set_config.get_python packages`

		- :func:`BacterialTyper.config.set_config.check_install_module`

		- :func:`BacterialTyper.config.extern_progs.min_python_module_version`

		- :func:`BacterialTyper.config.install_dependencies.python_package_install`

		- :func:`BacterialTyper.scripts.functions.file2dictionary`

	"""
	## get python packages installed
	my_packages_installed = get_python_packages(Debug)

	## debug messages
	if (Debug):
		print ("my_packages_installed :: ")
		print (my_packages_installed)

	## min versions for packages
	my_packages_requirements = extern_progs.min_python_module_version()

	## debug messages
	if (Debug):
		print ("my_packages_requirements")
		print (my_packages_requirements)

	## my module name conversion
	file_module_dependecies = extern_progs.file_list("module_dependencies")
	module_dependencies = functions.file2dictionary(file_module_dependecies, ',')

	## check each package
	for each in my_packages_requirements:
		## get min version
		min_version = my_packages_requirements[each]

		## get version installed in system
		installed = my_packages_installed[each] 

		## module name conversion
		module_name = module_dependencies[each]

		## debug messages
		if (Debug):
			print ("Module:", each)
			print ("Module name:", module_name)
			print ("Min_Version:", min_version)
			print (type(min_version))
			print ("Version installed:", installed)
			print (type(installed))

		## check if installed
		message = check_install_module(installed, module_name, min_version, 'Module')

		if (message == 'OK'):
			continue
		else:
			if (option_install):  # try to install
				if (Debug):
					print ("Install module: ", each)

				installed = install_dependencies.python_package_install(module_name, min_version)
				message2 = check_install_module(installed, module_name, min_version, 'Module')
				if (message2 == 'OK'):
					continue
				else:
					print ("+ Attent to install package: ", module_name, " failed. Install it manually to continue with BacterialTyper\n\n")
			else:
				print ("+ Please install manually package: ", module_name, " to continue with BacterialTyper\n\n")

################
def check_package_version(package, Debug):
	"""
	Retrieve python package version installed

	This is a modification of the original code from ARIBA (https://github.com/sanger-pathogens/ariba). 
	It basically uses pkg_resources.get_distribution(), pkg_resources.resource_filename() or imports module
	and retrieves version from __version__ variable.

	:param package: Python package name 
	:param Debug: True/False for debugging messages

	:type package: string
	:type Debug: boolean

	:returns: Version retrieved

	.. attention:: Be aware of Copyright

		The code implemented here was retrieved and modified from ARIBA (https://github.com/sanger-pathogens/ariba)

		Give them credit accordingly.
	"""

	try:
		version = pkg_resources.get_distribution(package).version
		if (Debug):
			print ("Method: pkg_resources.get_distribution(package).version")
	except:
		try:
			exec('import ' + package)
			version = eval(package + '.__version__')
			if (Debug):
				print ("Method: exec('import ' + package); version = eval(package + '.__version__')")
		except:
			try:
				if (Debug):
					print ("Method: pkg_resources.resource_filename(package, 'version.py')")
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
def get_perl_packages(Debug, file_name):
	"""
	Retrieves the version of the perl packages installed in the system.

	It retrieves the dependencies name conversion from file :file:`BacterialTyper/config/perl/perl_dependencies.csv`
	using function :func:`BacterialTyper.config.extern_progs.file_list` and :func:`BacterialTyper.scripts.functions.get_data`.
	For each module it retrieves the package version installed in the system using 
	:func:`BacterialTyper.config.set_config.check_perl_package_version`.	

	:returns: Dictionary containing for each perl module (key) the installed version (value).

	.. seealso:: This function relies on other ``BacterialTyper`` functions:

		- :func:`BacterialTyper.config.set_config.check_perl_package_version`

		- :func:`BacterialTyper.config.extern_progs.file_list` 

		- :func:`BacterialTyper.scripts.functions.get_data`

	"""
	## get info for perl modules
	perl_lib_dependecies_file = extern_progs.file_list(file_name)
	perl_lib_dependecies = functions.get_data(perl_lib_dependecies_file, ',', 'index_col=0')

	my_packages_installed = {}
	for index_name, row in perl_lib_dependecies.iterrows():
		module_name = row['module']
		installed = check_perl_package_version(module_name, Debug) ## check version installed in system
		if not (installed):
			installed = 'n.a.'
		my_packages_installed[index_name] = installed

	return (my_packages_installed)

##################
def check_perl_packages(file_name, Debug, option_install, install_path):
	"""
	Check the perl packages required

	This functions checks wether the packages installed in the system fulfilled the 
	minimum version specified in the configuration file. Details of the perl packages 
	required are available in :file:`BacterialTyper/config/perl/`. 

	It uses function :func:`BacterialTyper.config.set_config.get_perl_packages` to
	retrieve the version of the perl packages installed in the system. Then it uses
	:func:`BacterialTyper.config.extern_progs.min_perl_package_version` to retrieve the minimum
	version specified. It compares them using function :func:`BacterialTyper.config.set_config.check_install_module`.

	:param file_name: Name of the file to search within :file:`BacterialTyper/config/perl/`.
	:param Debug: True/False for debugging messages
	:param option_install: True/False for installing missing dependencies
	:param install_path: Install path for installing modules.

	:type file_name: string
	:type Debug: boolean
	:type option_install: boolean
	:type install_path: string

	:returns: Print messages if packages are installed.

	.. seealso:: This function relies on other ``BacterialTyper`` functions:

		- :func:`BacterialTyper.config.set_config.get_perl_packages`

		- :func:`BacterialTyper.config.set_config.check_install_module`

		- :func:`BacterialTyper.config.extern_progs.min_perl_package_version`

		- :func:`BacterialTyper.config.install_dependencies.perl_package_install`

	"""
	## get perl packages installed
	my_packages_installed = get_perl_packages(Debug, file_name)

	## debug messages
	if (Debug):
		print ("my_packages_installed :: ")
		print (my_packages_installed)

	## min versions for packages
	my_packages_requirements = extern_progs.min_perl_package_version(file_name)

	## debug messages
	if (Debug):
		print ("my_packages_requirements")
		print (my_packages_requirements)

	## get info for perl modules
	perl_lib_dependecies_file = extern_progs.file_list(file_name)
	perl_lib_dependecies = functions.get_data(perl_lib_dependecies_file, ',', 'index_col=0')

	## check each package
	for each in my_packages_requirements:
		## get min version
		min_version = my_packages_requirements[each]

		## get version installed in system
		installed = my_packages_installed[each] 

		## module name conversion
		module_name = perl_lib_dependecies.loc[each, 'module']

		## debug messages
		if (Debug):
			print ("Module:", each)
			print ("Module name:", module_name)
			print ("Min_Version:", min_version)
			print ("Version installed:", installed)
			
		## check if installed
		message = check_install_module(installed, module_name, min_version, 'Package')

		if (message == 'OK'):
			continue
		else:
			print (colored("** ATTENTION: Installation of perl modules is not supported",'red'))
			print ("+ Please install manually package: ", module_name, " to continue with BacterialTyper\n\n")

################
def check_perl_package_version(package, Debug):
	"""
	Retrieve perl package version installed

	It basically uses a one line perl command to load the package and print the version.

	:param package: package name 
	:param Debug: True/False for debugging messages

	:type package: string
	:type Debug: boolean

	:returns: Version retrieved
	"""

	perl_exe = get_exe('perl')
	perl_one_line_command = perl_exe + ' -M' + package + ' -e \'print $' + package + '::VERSION\';'

	if (Debug):
		print ("** DEBUG: perl command:\n")
		print (perl_one_line_command)

	## execute one line perl command
	output_one_line = functions.system_call(perl_one_line_command, True)
	return(functions.decode(output_one_line))

################
## IslandPath
################
def check_IslandPath(Debug, option_install, install_path):

	## get perl packages installed
	check_perl_packages("IslandPath_dependencies", Debug, option_install, install_path)

	## check additional software required
	print ("+ Check additional software for IslandPath optional analysis...")

################
## Miscellaneous
################
def print_module_comparison(module_name, message, color, tag):
	"""
	Creates print message for a given module, version, message

	:param module_name: Name of the module
	:param message: Message to include in the print message: OK | FAILED | NOT FOUND
	:param color: Print message color: green | orange | red
	:param tag: Tag to include: Module, package, software

	:type module_name: string
	:type message: string 
	:type color: string
	:type tag: string

	:returns: Print message
	"""
	print (colored("{:.<15}{:.>15}".format("%s: %s" %(module_name, tag), "[ %s ]" %message), color))

#########

##################
def check_install_module(installed, module_name, min_version, tag):
	"""
	Checks modules installation 

	Checks wether a module is installed and fulifilling requirements. 	
	It prints messages using :func:`BacterialTyper.config.set_config.print_module_comparison`.

	:param installed: Version string of the module installed.
	:param module_name: Module name
	:param min_version: Version string for the minimum version required

	:type installed: string
	:type module_name: string 
	:type min_version: string
	"""
	 ## Not installed
	if (installed == 'n.a.'):
		message = 'NOT FOUND'
		color = 'red'
		print_module_comparison(module_name, message, color, tag)

	# check version
	elif LooseVersion(installed) >= LooseVersion(min_version):
		message = 'OK'
		color = 'green'
		print_module_comparison(module_name, message, color, tag)

	else:
		message = 'FAILED'
		color = 'yellow'
		print_module_comparison(module_name, message, color, tag)

	## return message
	return (message)

#########
