#usr/bin/env python
'''
This module install external dependencies if not satistified
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
from BacterialTyper import extern_progs

##################
def install(software):
	print ("Install missing software: ", software)
	print ("To do....")
	
	## try to install: 
	print(colored("**Check paths or install it in the system and add it to $PATH environment variable.",'red'))

##################
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
