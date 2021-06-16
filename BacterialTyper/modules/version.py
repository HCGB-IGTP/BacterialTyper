#!/usr/bin/env python3
##########################################################
## Jose F. Sanchez										##
## Copyright (C) 2019 Lauro Sumoy Lab, IGTP, Spain		##
##########################################################
'''
Prints an index of version numbers for the different packages and other softwares employed here.
'''

## useful imports
from termcolor import colored
import os
import sys

## import my modules
from BacterialTyper.config import set_config
from BacterialTyper.config import extern_progs
import HCGB.functions.aesthetics_functions as HCGB_aes
from BacterialTyper import __version__ as pipeline_version

#########
def print_all():

	HCGB_aes.print_sepLine("+", 50, False)
	print ('Python:')
	HCGB_aes.print_sepLine("+", 50, False)
	print ('Python version:', str(sys.version))
	print ('\n')
	HCGB_aes.print_sepLine("+", 50, False)
	print ('Python packages:')
	extern_progs.print_package_version()
	HCGB_aes.print_sepLine("+", 50, False)
	print ('\n')

	HCGB_aes.print_sepLine("+", 50, False)
	print ('External dependencies:\n')
	HCGB_aes.print_sepLine("+", 50, False)
	extern_progs.print_dependencies()
	print ('\n')

	print ('Additional dependencies: databases, information, etc...')

	HCGB_aes.print_sepLine("*", 50, False)
	print ("ARIBA databases version..")
	HCGB_aes.print_sepLine("*", 50, False)
	print ("card -> ")
	print ("megares  ->")
	print ("plasmidfinder  ->")
	print ("resfinder  ->")
	print ("srst2_argannot  ->")
	print ("vfdb_core & vfdb_full  ->")
	print ("virulencefinder  ->")
	print ('\n')
	
	## kma databases
	## ...
	
#########
def only_us():
	print ("Version to be included when properly finished...\n\n")

	HCGB_aes.print_sepLine("+", 50, False)
	print ('BacterialTyper version: ' + pipeline_version)
	HCGB_aes.print_sepLine("+", 50, False)
	print ('\n')


#########
def run(options):
	HCGB_aes.pipeline_header("BacterialTyper", ver=pipeline_version)
	HCGB_aes.boxymcboxface("Version")

	if (options.option == 'all'):
		print (colored("\n+ BacterialTyper version:", 'yellow'))
		only_us()

		print (colored("\n+ Other softwares employed in the pipeline", 'yellow'))
		print_all()
		
	elif (options.option == 'only'):
		print (colored("\n+ BacterialTyper version:", 'yellow'))
		only_us()

	return()
