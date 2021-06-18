#!/usr/bin/env python3
##########################################################
## Jose F. Sanchez										##
## Copyright (C) 2019 Lauro Sumoy Lab, IGTP, Spain		##
##########################################################
"""
It generates a protein annotation using Prokka_ for the genomes assembled generated

Prokka_ contains several databases with known parameters for several kingdom and genus. Prokka option `--usegenus` is set as default.
Several options are available for kingdom and genus. See details below and see output from :code:`prokka --listdb`
shown in :func:`BacterialTyper.scripts.annotation.print_list_prokka`. If no genus matches user desired option, use option "Other".

.. include:: ../../links.inc	 	
"""

## useful imports
import time
import io
import os
import re
import sys
from sys import argv
from io import open
import pandas as pd
from termcolor import colored

## import my modules
import HCGB
import HCGB.functions.time_functions as HCGB_time
import HCGB.functions.system_call_functions as HCGB_sys
from BacterialTyper.config import set_config

#############################################
def print_list_prokka():
	"""
	Prints Prokka_ databases that has installed to use. It is the output from the call: 
	
	.. code-block:: sh

		prokka --listdb
	
	.. include:: ../../devel/results/print_list_prokka.txt
		:literal:
	
	.. seealso:: This function depends on other BacterialTyper functions called:
	
		- :func:`BacterialTyper.scripts.set_config.get_exe`
			
	.. include:: ../../links.inc	 	
	"""
	prokka_bin = set_config.get_exe('prokka')
	cmd = prokka_bin + " --listdb"
	HCGB_sys.system_call(cmd)

#############################################
def module_call(sequence_fasta, kingdom, genus, path, name, threads):
	"""
	Function that checks and generates annotation.
	
	- It uses Prokka_ via :func:`BacterialTyper.scripts.annotation.prokka_call`.
	
	- It checks if previously generated 
	
	- Once finished, it prints timestamp 
	
	:param sequence_fasta: Assembled sequences in fasta file format. 
	:param kingdom: Available kingdoms mode for Prokka software: Archaea|Bacteria|Mitochondria|Viruses
	:param genus: Available genus options for Prokka software. See details above.
	:param path: Absolute path to the output folder to include results.
	:param name: Sample name and tag to include in the annotation report and files.
	:param threads: Number of CPUs to use.
	  
	:type sequence_fasta: string
	:type kingdom: string
	:type genus: string 
	:type path: string 
	:type name: string 
	:type threads: integer 
	
	.. seealso:: This function depends on other BacterialTyper functions called:
	
		- :func:`BacterialTyper.scripts.set_config.get_exe`
		
		- :func:`HCGB.functions.time_functions.read_time_stamp`
		
		- :func:`HCGB.functions.time_functions.print_time_stamp`
				
		- :func:`HCGB.functions.time_functions.prokka_call`	

	.. include:: ../../links.inc	 	
	"""
	
	## check if previously assembled and succeeded
	filename_stamp = path + '/.success'

	if os.path.isdir(path):
		if os.path.isfile(filename_stamp):
			stamp =	HCGB_time.read_time_stamp(filename_stamp)
			print (colored("\tA previous command generated results on: %s [%s]" %(stamp, name), 'yellow'))
			return ()
	
	## call prokka
	prokka_bin = set_config.get_exe('prokka')
	dirname = prokka_call(prokka_bin, sequence_fasta, kingdom, genus, path, name, threads)

	## success stamps
	filename_stamp = path + '/.success'
	stamp =	HCGB_time.print_time_stamp(filename_stamp)

	return(dirname)	

#############################################
def prokka_call(prokka_bin, sequence_fasta, kingdom, genus, outdir_name, name, threads):
	"""Create system call for Prokka_ software. 
	
	It generates genome annotation using Prokka software. 
		
	:param prokka_bin: Path to the prokka binary file.
	:param sequence_fasta: Assembled sequences in fasta file format. 
	:param kingdom: Available kingdoms mode for Prokka software: Archaea|Bacteria|Mitochondria|Viruses
	:param genus: Available genus options for Prokka software. See details above.
	:param outdir_name: Absolute path to the output folder to include results.
	:param name: Sample name and tag to include in the annotation report and files.
	:param threads: Number of CPUs to use.
	  
	:type prokka_bin: string
	:type sequence_fasta: string
	:type kingdom: string
	:type genus: string 
	:type outdir_name: string 
	:type name: string 
	:type threads: integer 
	
	.. seealso:: Check description of output files generated in:
	
		- :ref:`Prokka-output-files`


	.. seealso:: This function depends on other BacterialTyper functions called:
	
		- :func:`BacterialTyper.scripts.set_config.get_exe`
		
		- :func:`HCGB.functions.time_functions.read_time_stamp`
		
		- :func:`HCGB.functions.time_functions.print_time_stamp`
				
		- :func:`BacterialTyper.scripts.annotation.prokka_call`	
	
	.. include:: ../../links.inc
	"""
	
	## set parameters and options for prokka
	print ("\n+ Starting annotation for: %s\n" %name)
	log_file = outdir_name + '/run.log'
	options = "--cdsrnaolap --addgenes --addmrna --kingdom " + kingdom
	if genus != "Other":
		options = options + " --usegenus --genus " + genus
	prokka = "%s --force --outdir %s --prefix %s --locustag %s %s --cpus %s %s 2> %s" %(prokka_bin, 
																					outdir_name, name, name, options, 
																					threads, sequence_fasta, log_file)
	HCGB_sys.system_call(prokka)
	return(outdir_name)

#############################################
def help_options():
	print ("\nUSAGE: python %s sequence_file path name CPUs kingdom prokka_bin\n"  %os.path.realpath(__file__))

#############################################
def main():

 	## control if options provided or help
	if len(sys.argv) > 1:
		print ("")
	else:
		help_options()
		exit()
	
	## argv
	seq_file = os.path.abspath(argv[1])
	folder = os.path.abspath(argv[2])
	name = argv[3]
	threads = int(argv[4])
	kingdom = argv[5]
	prokka_bin = argv[6]
	
	## call
	dir_annot=prokka_call(prokka_bin, seq_file, kingdom, folder, name, threads)
	print ("\n+ Annotation for sample %s has been generated in folder: %s" %(name, dir_annot))

######
'''******************************************'''
if __name__== "__main__":
	main()
