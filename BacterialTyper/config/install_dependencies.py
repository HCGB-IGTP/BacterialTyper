#!/usr/bin/env python3
##########################################################
## Jose F. Sanchez					##
## Copyright (C) 2019-2020 Lauro Sumoy Lab, IGTP, Spain	##
##########################################################
"""
Installs external dependencies if not satistified
"""
## this modules is an idea from ARIBA (https://github.com/sanger-pathogens/ariba)
## give credit to them appropiately

## [TODO]

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
from BacterialTyper.scripts import functions
from BacterialTyper.config import set_config
from BacterialTyper.config import extern_progs

##################
## [TODO]
## def install(software): edirect NCBI


## [TODO]
##################
def install(software, min_version, install_path):
	
	## install busco
	if (software == 'busco' or software == 'busco_plot'):
		install_BUSCO()
		
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
	print ("\nSome error ocurred while installing module [ %s ]." %module_name)
	print ("Path: " + path)
	print ("Please retry or install it mannually to continue with BacterialTyper")


##################
def perl_package_install(name, version2install, Debug):
	"""
	Installs perl package using CPAN
	
	This function installs packages archived in CPAN using cpan perl module.
	
	:param name: Perl package name
	:param version2install: Version to install
	:param Debug: True/False for debugging messages

	:type name: string
	:type version2install: string 
	:type Debug: boolean

	.. seealso: This function depends on other ``BacterialTyper`` functions such as:

		- :func:`BacterialTyper.config.install_dependencies.install_cpan_package`
	"""	
	print ()
	
##################	
def perl_package_install_targz(name, version2install, http_tar_gz, install_dir, Debug):
	"""
	Retrieves information for perl package installation
	
	This function installs packages archived in tar.gz files. It downloads and extracts
	``tar.gz`` files and install them locally. This option might not be the best option
	as long as it would not install all sub-dependencies for each package. We would use
	preferentially :func:`BacterialTyper.config.install_dependencies.perl_package_install_cpan`.  

	:param name: Perl package name
	:param version2install: Version to install
	:param http_tar_gz: FTP/https site of the tar gz perl package (cpan)
	:param install_dir: Installation directory
	:param Debug: True/False for debugging messages

	:type name: string
	:type version2install: string 
	:type http_tar_gz: string 
	:type install_dir: string 
	:type Debug: boolean

	.. seealso: This function depends on other ``BacterialTyper`` functions such as:

		- :func:`BacterialTyper.scripts.functions.wget_download`

		- :func:`BacterialTyper.scripts.functions.extract`

		- :func:`BacterialTyper.config.install_dependencies.install_package`

	"""	

	print (colored("Install missing perl package: " + name, 'yellow'))

	## create folders
	perlPackages = functions.create_subfolder("perl_packages", install_dir)
	path2download_name = functions.create_subfolder(name, perlPackages)

	## debugging messages
	if (Debug):
		print ("perlPackages: " + perlPackages)
		print ("path2download_name: " + path2download_name)

        ## download
	functions.wget_download(http_tar_gz, path2download_name)

	## get targz file
	tar_gz = functions.retrieve_matching_files(path2download_name, 'tar.gz')
	functions.extract(tar_gz[0], path2download_name)

	## debugging messages
	if (Debug):
		print ("** DEBUG: **")
		print ("http: " + http_tar_gz)
		print ('tar_gz: ')
		print (tar_gz)
		print ("*******************")

	## install
	path2download_extract_folder = os.path.join(path2download_name, name)
	blib_lib = install_package_targz(path2download_extract_folder, install_dir, Debug, name)

	## include in $PERL5LIB system variable
	print (colored("** ATTENTION:", 'yellow'))
	print ("Include the following path within your PERL5LIB variable:")
	print (blib_lib)
	
	## debugging messages
	if (Debug):
		print ("** DEBUG: **")
		print ("blib_lib: " + blib_lib)

	return(version2install)

##################

#######################
def install_package_targz(package_path, install_path, Debug, name):
	"""
	Install perl package downloaded
	
	:param package_path: Path where the perl package is extracted
	:param install_path: Path to install package

	:type package_path: string
	:type install_path: string

	:returns: Path where the module is installed to include in $PERL5LIB
	"""
	## change dir to package path
	os.chdir(package_path)

	print ("## Installing module: " + name + " ##")

	## perl Makefile.PL
	makefile_perl = functions.retrieve_matching_files(package_path, "Makefile.PL")
	perl_exe = set_config.get_exe("perl", Debug)

	## two installation possibilities: Makefile.PL or Build.PL
	if (makefile_perl):

		print ("+ Create make file")
		perl_MakeFile_cmd = perl_exe + ' ' + makefile_perl[0]

		## debug messages
		if (Debug):
			print ("** Debug: chdir " + package_path)
			print ("** Debug: perl_exe" + perl_exe)
			print ("** Debug: makefile_perl: " + makefile_perl[0])

		code_perl_make = functions.system_call(perl_MakeFile_cmd)
		code_make = ""

		##
		if (code_perl_make == 'OK'):
			## debug messages
			if (Debug):
				print ("** Debug: perl Makefile.PL successful")

			make_bin = set_config.get_exe("make", Debug)
			print ("+ Execute make file")
			code_make = functions.system_call(make_bin)

			if (code_make == 'OK'):
				if (Debug):
					print ("** Debug: make successful")
			else:
				print_error_message(name, package_path)
				return('n.a.')
		else:
			print_error_message(name, package_path)
			return('n.a.')
	else:
		## perl Makefile.PL
		Build_perl = functions.retrieve_matching_files(package_path, "Build.PL")
		if (Build_perl):
			print ("+ Create Build file")
			perl_build_cmd = perl_exe + ' ' + Build_perl[0]

			## debug messages
			if (Debug):
				print ("** Debug: chdir " + package_path)
				print ("** Debug: perl_exe" + perl_exe)
				print ("** Debug: Build_perl: " + Build_perl[0])

			code_perl_build = functions.system_call(perl_build_cmd)
			code_build = ""

			if (code_perl_build == 'OK'):
				## debug messages
				if (Debug):
					print ("** Debug: perl Build.PL successful")

				print ("+ Execute Build file")
				code_build = functions.system_call("./Build")

				if (code_build == 'OK'):
					if (Debug):
						print ("** Debug: build successful")
				else:
					print_error_message(name, package_path)
					return('n.a.')
			else:
				print_error_message(name, package_path)
				return('n.a.')
		else:
			print_error_message(name, package_path)
			return('n.a')


	## copy files an finish installation
	print ("+ Include files into the PERL path")
	blib_lib = os.path.join(package_path, 'blib', 'lib')
	if os.path.isdir(blib_lib):
		return(blib_lib)
	else:
		print_error_message(name, package_path)
		return ('n.a.')


def python_package_install(package, version2install):
	print (colored("Install missing python package: " + package, 'yellow'))
	versioninstalled='0.1'
	return (versioninstalled)

##################
def install_BUSCO():
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


	install_path = os.path.abspath(sys.argv[1])
	module_path = os.path.abspath(sys.argv[2])
	name = sys.argv[3]

	install_package(module_path, install_path, False, name)


if __name__ == "__main__":
	main()
