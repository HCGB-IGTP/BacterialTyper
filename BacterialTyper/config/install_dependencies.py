#!/usr/bin/env python3
##########################################################
## Jose F. Sanchez					##
## Copyright (C) 2019-2020 Lauro Sumoy Lab, IGTP, Spain	##
##########################################################
"""
Installs external dependencies if not satisfied
"""
## this modules is an idea from ARIBA (https://github.com/sanger-pathogens/ariba)
## give credit to them appropriately

## [TODO]

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

## import my modules
from BacterialTyper.scripts import functions
from BacterialTyper.config import set_config
from BacterialTyper.config import extern_progs

##################
## [TODO]
## def install(software): edirect NCBI


## [TODO]
##################
def get_info_software():
	"""Read software information
	
	Reads information stored in file :file:`BacterialTyper.config.software.software_details.csv`
	and returns pandas data frame.	
	"""
	info_file = os.path.abspath(os.path.join(os.path.dirname( __file__ ), 'software', 'software_details.csv'))
	return(functions.get_data(info_file, ',', 'index_col=0'))

##################
def install(software, min_version, install_path):
	
	##
	Debug=True
	
	## install busco
	if (software == 'busco' or software == 'busco_plot'):
		install_BUSCO(install_path)
	
	## install edirect
	elif (software == 'efetch' or software == 'esearch' or software == 'xtract' ):
		install_edirect(install_path)
	else:
		install_soft(software, min_version, install_path, Debug)
			
	print ("Install missing software: ", software)
	print ("To do....")
	
	## try to install: 
	print(colored("**Check paths or install it in the system and add it to $PATH environment variable.",'green'))

	versionInstalled = 'n.a.'
	return (versionInstalled)

#######################
def print_error_message(module_name, path):
	"""
	Print error messages generated during installation

	:param module_name: 
	:param path:

	:type module_name: string 
	:type path: string

	:returns: Print messages.
	"""
	print ("\nSome error happened while installing module [ %s ]." %module_name)
	print ("Path: " + path)
	print ("Please retry or install it manually to continue with BacterialTyper")

#############################3
def python_package_install(package, version2install):
	print (colored("Install missing python package: " + package, 'yellow'))
	versioninstalled='0.1'
	return (versioninstalled)

##################

#######################
def install_soft(software, min_version, install_path, Debug):
	
	# check info from file: software_details.txt
	info = get_info_software()

	## print debugging messages
	if Debug:
		print(colored("** Debug: info data frame retrieved from file", 'yellow'))
		print(info)	
	
	## get info for file: web site & ext
	type= info.loc[software, 'type']

	if (type=='extract'):
		install_binary(software, min_version, install_path, Debug)
	elif (type=='git'):
		install_git_repo(software, min_version, install_path, Debug)
	else:
		print()

#######################
def install_git_repo(software, min_version, install_path, Debug):
	""" """
	# check info from file: software_details.txt
	info = get_info_software()

	## get info for file: web site & ext
	git_repo= info.loc[software, 'site']
	folder_name = info.loc[software, 'folder']
	bin_name = info.loc[software, 'bin_name']
	
	print ("+ Installing:")
	print ("\tSoftware: ", software)
	print ("\tPath: ", install_path)
	print ("\tGit repo: ", git_repo)
	
	## debug messages
	if Debug:
		print(colored("** Debug:", 'yellow'))
		print(colored("\nGit repo: %s" %git_repo, 'yellow'))
		print(colored("Folder clone: %s" %folder_name, 'yellow'))
		print(colored("Binary name: %s\n" %bin_name, 'yellow'))
	
	## option
	option = info.loc[software, 'ext']
	folder_name = info.loc[software, 'folder']
	binary_name = info.loc[software, 'bin_name']
	folder_path = os.path.join(install_path, folder_name)

	## current path
	current_path = os.getcwd()
	
	## git clone repo
	print ('+ Clone repository...')
	git_exe = set_config.get_exe('git', Debug)
	
	if os.path.exists(folder_path):
		## pull
		os.chdir(folder_path)
		cmd = git_exe + ' pull'
	else:
		## clone
		cmd = git_exe + ' clone ' + git_repo 
	
	## call git
	functions.system_call(cmd)

	## Compile
	print ('+ Compile software...')
	if (option == 'make'):
		## make
		os.chdir(folder_path)
		functions.system_call('make')
		
	else:
		## no need to compile
		print ()
	
	## get software path
	file_software = os.path.join(folder_path, binary_name)
	
	## add some exceptions
	if (software == 'prokka'):
		## get software path
		file_software = os.path.join(folder_path, 'bin', binary_name)
		
	## debug messages
	if Debug:
		print(colored("** Debug:", 'yellow'))
		print(colored('folder_path: %s' %folder_path, 'yellow'))
		print(colored('file_software: %s' %file_software, 'yellow'))
	
	## change permissions
	functions.chmod_rights(file_software, stat.S_IRWXU) ## execute permissions for group
	
	## get version
	VersionInstalled = set_config.get_version(software, file_software, Debug)
	
	## debug messages
	if Debug:
		print(colored("** Debug: VersionInstalled: %s" %VersionInstalled, 'yellow'))
	
	## chdir to previous path
	os.chdir(current_path)
	
	##
	return (VersionInstalled)
	
	
#######################
def install_binary(software, min_version, install_path, Debug):
	"""Install binary software
	
	For some software packages there are already compiled binaries for Linux. 
	This function retrieves them from the available URL, extracts them and 
	sets the appropriate system path for further usage. 
	
	.. seealso:: This function depends on other ``BacterialTyper`` functions:

		- :func:`BacterialTyper.config.install_dependencies.get_info_software`
		
		- :func:`BacterialTyper.config.set_config.get_version`
		
		- :func:`BacterialTyper.config.set_config.get_version`
		
	"""
	
	# check info from file: software_details.txt
	info = get_info_software()

	## get info for file: web site & ext
	site= info.loc[software, 'site']
	ext = info.loc[software, 'ext']
	folder_name = info.loc[software, 'folder']
	bin_name = info.loc[software, 'bin_name']
	version2install = info.loc[software, 'version']
	
	print ("+ Installing:")
	print ("\tSoftware: ", software)
	print ("\tVersion: ", version2install)
	print ("\tPath: ", install_path)
	
	## debug messages
	if Debug:
		print(colored("\nSite: %s" %site, 'yellow'))
		print(colored("Extension: %s" %ext, 'yellow'))
		print(colored("Folder downloaded: %s" %folder_name, 'yellow'))
		print(colored("Binary name: %s\n" %bin_name, 'yellow'))
	
	## check if folder already available
	folder_software = os.path.join(install_path, folder_name)
	if os.path.exists(folder_software):
		shutil.rmtree(folder_software)
			
	## check if file already available
	compress_file_name = os.path.join(install_path, site.rsplit('/', 1)[-1])
	if not functions.is_non_zero_file(compress_file_name):
		## download
		functions.wget_download(site, install_path)
	
	## extract
	functions.extract(compress_file_name, install_path, remove=False)
	file_software = os.path.join(folder_software, bin_name)
	
	## debug messages
	if Debug:
		print(colored("** Debug:", 'yellow'))
		print(colored('file2extract: %s' %compress_file_name, 'yellow'))
		print(colored('folderExtracted: %s' %folder_software, 'yellow'))
		print(colored('file_software: %s' %file_software, 'yellow'))
	
	## change permissions
	functions.chmod_rights(file_software, stat.S_IRWXU) ## execute permissions for group
	
	## get version
	VersionInstalled = set_config.get_version(software, file_software, Debug)
	
	## debug messages
	if Debug:
		print(colored("** Debug: VersionInstalled: %s" %VersionInstalled, 'yellow'))

	## return	
	return(VersionInstalled)
	
##################
def install_edirect(install_path):
	"""
	Installs and configures Edirect 
	
	Read further information of the Edirect utilities in https://www.ncbi.nlm.nih.gov/books/NBK179288/
	"""
	
	print ()
	
##################
def install_BUSCO(install_path):
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


def main():


	install_soft('spades', '', '/home/jfsanchez/software', True)
	install_soft('prokka', '', '/home/jfsanchez/software', True)
	install_soft('trimmomatic', '', '/home/jfsanchez/software', True)

if __name__ == "__main__":
	main()


