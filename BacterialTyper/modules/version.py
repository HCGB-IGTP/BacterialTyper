#!/usr/bin/env python3
'''
This code generates prints an index of version numbers for the different packages and other softwares employed here.
Jose F. Sanchez
Copyright (C) 2019 Lauro Sumoy Lab, IGTP, Spain
'''

## import my modules
from BacterialTyper import config, functions
from termcolor import colored
import os, sys

from BacterialTyper import extern_progs
#from BacterialTyper import __version__ as ariba_version ## to include when distributed
pipeline_version = "0.0.2"

## ToDo set automatic from pip list
python_packages_BacterialTyper = ('ariba', 'bs4', 'dendropy', 'pyfastaq', 'pymummer', 'pysam')

#########
def print_all():

	functions.print_sepLine("+", 50)
	print ('BacterialTyper version: ' + pipeline_version)
	functions.print_sepLine("+", 50)
	print ('\n')

	functions.print_sepLine("+", 50)
	print ('External dependencies:\n')
	functions.print_sepLine("+", 50)
	extern_progs.dependencies()
	print ('\n')

	functions.print_sepLine("+", 50)
	print ('Python:')
	functions.print_sepLine("+", 50)
	print ('Python version:', str(sys.version))
	print ('\n')
	functions.print_sepLine("+", 50)
	print ('Python packages:')
	check_python_packages(python_packages_BacterialTyper)        
	functions.print_sepLine("+", 50)
	print ('\n')

	print ('Additional dependencies: databases, information, etc...')

	functions.print_sepLine("*", 50)
	print ("ARIBA databases version..")
	functions.print_sepLine("*", 50)
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

#########
def run(options):
	if (options.option == 'all'):
		print (colored("\n+ BacterialTyper citation:", 'blue'))
		only_us()

		print (colored("\n+ Other softwares employed in the pipeline", 'blue'))
		print_all()
		
	elif (options.option == 'only'):
		print (colored("\n+ BacterialTyper version:", 'blue'))
		only_us()

#########
def check_python_packages(givenList):
	## this function is from ARIBA (https://github.com/sanger-pathogens/ariba)
	## give credit to them appropiately

	for package in givenList:
		try:
			exec('import ' + package)
			version = eval(package + '.__version__')
			path = eval(package + '.__file__')
			print (package + ' version: ' + version)
    
		except:
			version = 'NOT_FOUND'
			path = 'NOT_FOUND'
			print (package + ' version: ' + version)
	
