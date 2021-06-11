#!/usr/bin/env python3
##########################################################
## Jose F. Sanchez					##
## Copyright (C) 2019-2020 Lauro Sumoy Lab, IGTP, Spain	##
##########################################################
"""
Installs external dependencies if not satisfied
"""

## useful imports
import os
import io
import sys
import re
import shutil
import urllib
from io import open
from sys import argv
import stat
import subprocess
import pandas as pd
from termcolor import colored
from distutils.version import LooseVersion

## import my modules
from BacterialTyper.config import set_config
from BacterialTyper.config import extern_progs

import HCGB.functions.aesthetics_functions as HCGB_aes
import HCGB.functions.time_functions as HCGB_time
import HCGB.functions.main_functions as HCGB_main
import HCGB.functions.files_functions as HCGB_files
import HCGB.functions.system_call_functions as HCGB_sys

##################
def get_info_software():
	"""Read software information
	
	Reads information stored in file :file:`BacterialTyper.config.software.software_details.csv`
	and returns pandas data frame.	
	"""
	info_file = os.path.abspath(os.path.join(os.path.dirname( __file__ ), 'software', 'software_details.csv'))
	return(HCGB_main.get_data(info_file, ',', 'index_col=0'))

#######################
def print_error_message(module_name, path, option):
	"""
	Print error messages generated during installation

	:param module_name: 
	:param path:

	:type module_name: string 
	:type path: string

	:returns: Print messages.
	"""
	print ("\nSome error happened while installing %s [ %s ]." %(option, module_name))
	print ("Path: " + path)
	print ("Please retry or install it manually to continue with BacterialTyper")

#############################3
def python_package_install(package, version2install):
	print (colored("Install missing python package: " + package, 'yellow'))
	versioninstalled='0.1'
	return (versioninstalled)

##################
def get_install_R_files():
	install_R = os.path.abspath(os.path.join(os.path.dirname( __file__ ), 'R', 'install_package.R'))
	install_github_package = os.path.abspath(os.path.join(os.path.dirname( __file__ ), 'R', 'install_github.R'))
	return(install_R, install_github_package)

##################
def install_R_packages(package, source, install_path, extra):
	
	(install_R, install_github_package) = get_install_R_files()
	
	HCGB_files.create_folder(install_path)
	Rscript_exe = set_config.get_exe('Rscript')
	print("+ Installing %s package..." %package)
	install_file = install_R
	if (source == 'github'):
		install_file = install_github_package
		package= extra + '/' + package
	
	cmd_R = '%s %s -l %s -p %s' %(Rscript_exe, install_file, package, install_path)
	HCGB_sys.system_call(cmd_R)
	
	## check if exists or try to install
	MLSTar_package = os.path.join(install_path, 'MLSTar')
	if os.path.exists(MLSTar_package):
		RDir_package = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'R', 'R_package.info.txt')
		HCGB_main.printList2file(RDir_package, [install_path])
	else:
		print_error_message(package, "No R package found", 'package')
		print ('Please install manually to proceed...')

##################
def install(software, min_version, install_path, Debug):
	
	(path2Export, versionInstalled) = install_soft(software, min_version, install_path, Debug)
			
	## failed to install:
	if not path2Export:
		print(colored("**Check paths or install it in the system and add it to $PATH environment variable.",'yellow'))
		return ()
	
	else:
		## add to $PATH: include in environment bin
		env_bin_directory = os.path.dirname(os.environ['_'])
		
		print ("\n+ Add software to path")

		file_list = []

		## unique file to export
		if (software == 'fastqc' or software == 'trimmomatic'):
			file_list.append(path2Export)
		
		else:
		## all folder
			if (software == 'spades'):
				pathToExport = os.path.join(path2Export, 'bin')
		
			if (software == 'prokka'):
				pathToExport = os.path.join(path2Export, 'bin')
			
			file_list = HCGB_main.get_fullpath_list(path2Export)
			
			## add binaries compiled for linux
			if (software == 'prokka'):
				pathToExport2 = os.path.join(path2Export, 'binaries', 'linux')
				file_list = file_list + HCGB_main.get_fullpath_list(pathToExport2)
			
		## discard some files obtain
		file_list = [s for s in file_list if '.a' not in s]
		file_list = [s for s in file_list if '.c' not in s]
		file_list = [s for s in file_list if '.o' not in s]
		file_list = [s for s in file_list if '.h' not in s]		
		file_list = [s for s in file_list if '.git' not in s]
		file_list = [s for s in file_list if '.git/' not in s]
		file_list = [s for s in file_list if '.gitignore' not in s]
		file_list = [s for s in file_list if 'Makefile' not in s]
		file_list = [s for s in file_list if '.pdf' not in s]
		file_list = [s for s in file_list if '.tar.gz' not in s]
		file_list = [s for s in file_list if 'README.md' not in s]
		file_list = [s for s in file_list if '__pycache__' not in s]
		file_list = [s for s in file_list if 'db/' not in s]
		file_list = [s for s in file_list if 'doc/' not in s]
		file_list = [s for s in file_list if 'test/' not in s]
		file_list = [s for s in file_list if 'aux/' not in s]
			
		## debug messages
		if Debug:
			print(colored("** Debug: list to include in path",'yellow'))
			print (file_list)
			print()
		
		## create symbolic link in bin directory in environment
		HCGB_main.get_symbolic_link(file_list, env_bin_directory)
		print(colored("**Software (%s - Version: %s) installed in the system and add it to $PATH environment variable." %(software, versionInstalled),'green'))

	return (versionInstalled)

#######################
def install_soft(software, min_version, install_path, Debug):
	# check info from file: software_details.txt
	info = get_info_software()
		
	## install busco
	if (software == 'busco' or software == 'busco_plot'):
		code = install_BUSCO(install_path)
		folder_software=""
		software='busco'
	else:		
		## install edirect
		if (software == 'efetch' or software == 'esearch' or software == 'xtract' ):
			software = 'edirect'			
		## install blast
		elif (software == 'tblastn' or software == 'blastn' or software == 'makeblastdb' ):
			software = 'blast'
		
		## check
		if software not in info.index:
			print_error_message(software, "No dependencies found", 'software')
			return ('','')
		
		## get info for file: web site & ext
		site= info.loc[software, 'site']
		ext = info.loc[software, 'ext']
		folder_name = info.loc[software, 'folder']
		bin_name = info.loc[software, 'bin_name']
		version2install = info.loc[software, 'version']
		type= info.loc[software, 'type']
			
		## check if folder already available
		folder_software = os.path.join(install_path, folder_name)
	
		## print debugging messages
		if Debug:
			print(colored("** Debug: info data frame retrieved from file", 'yellow'))
			print(info)	
		
		## Git clone, extract or compile according to software 
		if (type=='extract'):
			print ("+ Extracting binary:")
			print ("\tSoftware: ", software)
			print ("\tURL: ", site)
			print ("\tVersion: ", version2install)
			print ("\tPath: ", install_path)

			code = install_binary(folder_software, site, install_path, Debug)
		
		elif (type=='git'):
			print ("+ Installing from git:")
			print ("\tSoftware: ", software)
			print ("\tPath: ", install_path)
			print ("\tGit repo: ", site)
			
			code = install_git_repo(site, folder_software, install_path, ext, Debug)	
		
		elif (type=='source'):
			print ("+ Installing from source:")
			print ("\tSoftware: ", software)
			print ("\tURL: ", site)
			print ("\tVersion: ", version2install)
			print ("\tPath: ", install_path)
			
			code = install_source(software, install_path, Debug)
		else:
			print()
	
	## some error ocurred
	if not code:
		print_error_message(software, folder_software, 'software')
		return ('','')
	
	## check installation		
	VersionInstalled = ""
	
	## add exceptions: blast
	if (software == 'blast'):
		## get version as an example using blastn
		folder_software = os.path.join(folder_software, bin_name)
		file_software = os.path.join(folder_software, 'blastn')
		HCGB_sys.chmod_rights(file_software, stat.S_IRWXU) ## execute permissions for group
		VersionInstalled = set_config.get_version('blastn', file_software, Debug)
		
	elif (software == 'edirect'):
		## run setup.sh
		file_setup = os.path.join(folder_software, 'setup.sh')
		HCGB_sys.system_call(file_setup)
		
		## get version as an example using efetch
		file_software = os.path.join(folder_software, 'efetch')
		HCGB_sys.chmod_rights(file_software, stat.S_IRWXU) ## execute permissions for group
		VersionInstalled = set_config.get_version('efetch', file_software, Debug)
	
	elif (software == 'busco'):
		## setup busco
		print()

	else:
		## get version
		file_software = os.path.join(folder_software, bin_name)
		HCGB_sys.chmod_rights(file_software, stat.S_IRWXU) ## execute permissions for group
		VersionInstalled = set_config.get_version(software, file_software, Debug)

	## debug messages
	if Debug:
		print(colored("** Debug:", 'yellow'))
		print(colored('folderExtracted: %s' %folder_software, 'yellow'))
		print(colored('file_software: %s' %file_software, 'yellow'))
		print(colored("VersionInstalled: %s" %VersionInstalled, 'yellow'))
	
	if (VersionInstalled == 'na'):
		print_error_message(software, folder_software, 'software')
		return ('','')
	
	## check version
	if min_version:
		## debug messages
		if Debug:
			print(colored("** Debug: Min version required: %s" %min_version, 'yellow'))
	
		## check min version
		if LooseVersion(VersionInstalled) >= LooseVersion(min_version):
			## debug messages
			if Debug:
				print(colored("** Debug: Min version satisfied: %s > %s" %(VersionInstalled, min_version), 'yellow'))
		else:
			## min version not satisfied
			## debug messages
			if Debug:
				print(colored("** Debug: Min version not satisfied: %s < %s" %(VersionInstalled, min_version), 'red'))
			
			## error
			return ('','')
	
	## return path to include in $PATH
	if (software == 'fastqc' or software == 'trimmomatic'):
		pathToExport = file_software
	
	elif (software == 'spades'):
		pathToExport = os.path.join(folder_software, 'bin')
	else:
		pathToExport = folder_software

	if (software == 'busco'):
		print()

	if (software == 'prokka'):
		## it is necessary to setup the database
		setup_cmd = file_software + ' --setupdb'
		print ("+ Set up prokka databases...")
		HCGB_sys.system_call(setup_cmd)
			## return
	
	return(pathToExport, VersionInstalled)

#######################
def install_git_repo(git_repo, folder_sofware, install_path, option, Debug):
	""" """
	
	## current path
	current_path = os.getcwd()
	os.chdir(install_path)
	
	## git clone repo
	print ('+ Using git to get code...')
	git_exe = set_config.get_exe('git', Debug)
	
	if os.path.exists(folder_sofware):
		print ('+ Clone repository...')
		## pull
		os.chdir(folder_sofware)
		cmd = git_exe + ' pull'
	else:
		print ('+ Clone repository...')
		## clone
		cmd = git_exe + ' clone ' + git_repo 
	
	## call git
	HCGB_sys.system_call(cmd)

	## compile if necessary
	if (option == 'make'):
		## Compile
		print ('+ Compile software...')
		## make
		os.chdir(folder_sofware)
		HCGB_sys.system_call('make')
	
	## chdir to previous path
	os.chdir(current_path)
	
	return(True)	

#######################
def install_binary(folder_software, site, install_path, Debug):
	"""Install binary software
	
	For some software packages there are already compiled binaries for Linux. 
	This function retrieves them from the available URL, extracts them and 
	sets the appropriate system path for further usage. 
	
	.. seealso:: This function depends on other ``BacterialTyper`` functions:

		- :func:`BacterialTyper.config.install_dependencies.get_info_software`
		
		- :func:`BacterialTyper.config.set_config.get_version`
		
		- :func:`BacterialTyper.config.set_config.get_version`
		
	"""

	if os.path.exists(folder_software):
		shutil.rmtree(folder_software)
			
	## check if file already available
	compress_file_name = os.path.join(install_path, site.rsplit('/', 1)[-1])
	if not HCGB_main.is_non_zero_file(compress_file_name):
		## download
		HCGB_sys.wget_download(site, install_path)
	
	## extract
	HCGB_files.extract(compress_file_name, install_path, remove=False)
	return(True)

##################
def install_source(software, install_path, Debug):
	return()
	
	return(True)

##################
def install_BUSCO(install_path):
	
	return(False)
	## git clone https://gitlab.com/ezlab/busco.git

	## install within our environment
	## run python setup.py install --user 
	## export PYTHONPATH=$PYTHONPATH:/path/to/busco/build/lib

	## add BUSCO_CONFIG_FILE to path
	folder = ""
	fileConfig = BUSCO_config(folder)	
	#export BUSCO_CONFIG_FILE=fileConfig"

	# AUGUSTUS_CONFIG_PATH
	#cp -r /path/to/AUGUSTUS/augustus-3.2.3/config /folder/augustus/config
	#export AUGUSTUS_CONFIG_PATH="/folder/augustus/config"
	
	return(True)

##################
def BUSCO_config():

	## set BUSCO_CONFIG_FILE on $PATH
	folder = ""
	config_file = folder + '/config_BUSCO.ini'
	file_hd = open(config_file, 'w')

	## write BUSCO configuration into file	
	file_hd.write("# BUSCO specific configuration")
	file_hd.write("# It overrides default values in code and dataset cfg, and is overridden by arguments in command line")
	file_hd.write("# Uncomment lines when appropriate")
	file_hd.write("[busco]")
	file_hd.write("# Input file")
	file_hd.write(";in = ./sample_data/target.fa")
	file_hd.write("# Run name, used in output files and folder")
	file_hd.write(";out = SAMPLE")
	file_hd.write("# Where to store the output directory")
	file_hd.write(";out_path = ./sample_data")
	file_hd.write("# Path to the BUSCO dataset")
	file_hd.write(";lineage_path = ./sample_data/example")
	file_hd.write("# Which mode to run (genome / protein / transcriptome)")
	file_hd.write(";mode = genome")
	file_hd.write("# How many threads to use for multithreaded steps")
	file_hd.write(";cpu = 1")
	file_hd.write("# Domain for augustus retraining, eukaryota or prokaryota")
	file_hd.write("# Do not change this unless you know exactly why !!!")
	file_hd.write(";domain = eukaryota")
	file_hd.write("# Force rewrite if files already exist (True/False)")
	file_hd.write(";force = False")
	file_hd.write("# Restart mode (True/False)")
	file_hd.write(";restart = False")
	file_hd.write("# Blast e-value")
	file_hd.write(";evalue = 1e-3")
	file_hd.write("# Species to use with augustus, for old datasets only")
	file_hd.write(";species = fly")
	file_hd.write("# Augustus extra parameters")
	file_hd.write("# Use single quotes, like this: '--param1=1 --param2=2'")
	file_hd.write(";augustus_parameters = ''")
	file_hd.write("# Tmp folder")
	file_hd.write(";tmp_path = ./tmp/")
	file_hd.write("# How many candidate regions (contigs, scaffolds) to consider for each BUSCO")
	file_hd.write(";limit = 3")
	file_hd.write("# Augustus long mode for retraining (True/False)")
	file_hd.write(";long = False")
	file_hd.write("# Quiet mode (True/False)")
	file_hd.write(";quiet = False")
	file_hd.write("# Debug logs (True/False), it needs Quiet to be False")
	file_hd.write(";debug = True")
	file_hd.write("# tar gzip output files (True/False)")
	file_hd.write(";gzip = False")
	file_hd.write("# Force single core for the tblastn step")
	file_hd.write(";blast_single_core = True\n")
	#
	file_hd.write("# You will need to set the BUSCO_CONFIG_FILE environment variable to define a custom path ")
	file_hd.write("# (including the filename) to your config.ini file.")
	file_hd.write("# e.g. ")
	file_hd.write("# export BUSCO_CONFIG_FILE=~/busco/config.ini")
	file_hd.write("#")
	file_hd.write("# For full details of how to configure the config.ini consult the userguide in the directory above this")
	file_hd.write("# and the example config.ini in thei directory.")
	file_hd.write("#")
	file_hd.write("# The following paths should be correct\n")
	#
	file_hd.write("[tblastn]")
	file_hd.write("# path to tblastn")
	file_hd.write("path = /software/debian-8/bin/\n")
	#
	file_hd.write("[makeblastdb]")
	file_hd.write("# path to makeblastdb")
	file_hd.write("path = /software/debian-8/bin/\n")
	#
	file_hd.write("[augustus]")
	file_hd.write("# path to augustus")
	file_hd.write("path = /soft/bio/augustus/bin/\n")
	#
	file_hd.write("[etraining]")
	file_hd.write("# path to augustus etraining")
	file_hd.write("path = /soft/bio/augustus/bin/\n")
	#
	file_hd.write("# path to augustus perl scripts, redeclare it for each new script\n")
	file_hd.write("[gff2gbSmallDNA.pl]")
	file_hd.write("path = /soft/bio/augustus/scripts/")
	file_hd.write("[new_species.pl]")
	file_hd.write("path = /soft/bio/augustus/scripts/")
	file_hd.write("[optimize_augustus.pl]")
	file_hd.write("path = /soft/bio/augustus/scripts/\n")
	#
	file_hd.write("[hmmsearch]")
	file_hd.write("# path to HMMsearch executable")
	file_hd.write("path = /software/debian-8/bin/\n")
	#
	file_hd.write("[Rscript]")
	file_hd.write("# path to Rscript, if you wish to use the plot tool")
	file_hd.write("path = /software/debian-8/bin/\n")
	#
	file_hd.write("# You will also need to set a variable to your own augustus config directory using the variable AUGUSTUS_CONFIG_PATH")
	file_hd.write("# eg export AUGUSTUS_CONFIG_PATH=~/augustus")
	file_hd.write("# You can find example config files in the directory /soft/bio/augustus/config.\n")
	#
	file_hd.close()

##################
def main():
	""" Just for debugging purposes"""
	
	set_config.check_R_packages(False, True, '')
	exit()
	
	install('fastqc', '0.11.9', '/home/jfsanchez/software', True)
	install('spades', '3.14', '/home/jfsanchez/software', True)
	install('prokka', '', '/home/jfsanchez/software', True)
	install('trimmomatic', '0.39', '/home/jfsanchez/software', True)
	install('kma', '', '/home/jfsanchez/software', True)
	install('blastn', '2.10', '/home/jfsanchez/software', True)
	install('efetch', '11.7', '/home/jfsanchez/software', True)

##################
if __name__ == "__main__":
	main()


