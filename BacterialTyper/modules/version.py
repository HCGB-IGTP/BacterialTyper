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
from BacterialTyper.scripts import functions


#from BacterialTyper.scripts import __version__ as ariba_version ## to include when distributed
pipeline_version = "0.0.2"

#########
def print_all():

	functions.print_sepLine("+", 50, False)
	print ('Python:')
	functions.print_sepLine("+", 50, False)
	print ('Python version:', str(sys.version))
	print ('\n')
	functions.print_sepLine("+", 50, False)
	print ('Python packages:')
	extern_progs.print_package_version()
	functions.print_sepLine("+", 50, False)
	print ('\n')

	functions.print_sepLine("+", 50, False)
	print ('External dependencies:\n')
	functions.print_sepLine("+", 50, False)
	extern_progs.print_dependencies()
	print ('\n')

	print ('Additional dependencies: databases, information, etc...')

	functions.print_sepLine("*", 50, False)
	print ("ARIBA databases version..")
	functions.print_sepLine("*", 50, False)
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

	functions.print_sepLine("+", 50, False)
	print ('BacterialTyper version: ' + pipeline_version)
	functions.print_sepLine("+", 50, False)
	print ('\n')


#########
def run(options):
	functions.pipeline_header()
	functions.boxymcboxface("Version")

	if (options.option == 'all'):
		print (colored("\n+ BacterialTyper version:", 'yellow'))
		only_us()

		print (colored("\n+ Other softwares employed in the pipeline", 'yellow'))
		print_all()
		
	elif (options.option == 'only'):
		print (colored("\n+ BacterialTyper version:", 'yellow'))
		only_us()

	return()
