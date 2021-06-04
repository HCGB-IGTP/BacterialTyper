#!/usr/bin/env python3
##########################################################
## Jose F. Sanchez										##
## Copyright (C) 2019 Lauro Sumoy Lab, IGTP, Spain		##
##########################################################
"""
Calls multiQC_ to generate HTML statistics reports.

.. include:: ../../links.inc 
"""
## useful imports
import os
import io
import sys
from io import open
from sys import argv
from termcolor import colored

## import my modules
import HCGB
import HCGB.functions.system_call_functions as HCGB_sys
import HCGB.functions.main_functions as HCGB_main
from BacterialTyper.config import set_config

############
def multiQC_module_call(givenList, name, path, option):
	"""
	Prepares files for multiQC report generation.
	
	:param givenList: List of folder to search for multiQC report.
	:param name: Name to include in the html report.
	:param path: Absolute path for the output folder.
	:param option: Some options to provide to multiQC_call.
	
	:type givenList: list
	:type name: string
	:type path: string
	:type option: string
	
	.. seealso:: This function depends on other BacterialTyper functions called:
	
		- :func:`HCGB.functions.main_functions.printList2file`
		
		- :func:`BacterialTyper.scripts.multiQC_report.multiQC_call`
	
	"""
	pathFile = path + '/' + 'samples.txt'
	HCGB_main.printList2file(pathFile, givenList)
	multiQC_call(pathFile, name, path, option)	
	
############
def multiQC_call(pathFile, name, folder, option):
	"""
	multiQC_ report generation call.
	
	:param pathFile: File containing list of files to include in report.
	:param name: Name to include in the html report.
	:param folder: Absolute path for the output folder.
	:param option: Options to provide to multiQC call.
	
	:type pathFile: string
	:type name: string 
	:type folder: string 
	:type option: string
	
	:returns: :func:`BacterialTyper.scripts.functions.system_call` output (OK/FALSE)
		
	.. seealso:: This function depends on other BacterialTyper functions called:
	
		- :func:`HCGB.functions.system_call_functions`
	
	"""
	multiqc_bin = "multiqc" ## if we activate the environment it should be in $PATH
	## set options for call
	cmd = "%s --force -o %s -n %s -l %s -p -i 'MultiQC report' -b 'HTML report generated for multiple samples and steps' %s" %(multiqc_bin, folder, name, pathFile, option)
	
	## if a report was previously generated in the folder 
	## force to delete and generate a new one
	return(HCGB_sys.system_call(cmd))

############
def	help_options():
	"""
	multiQC_ help options for call as a single script.
	"""
	print ("\nUSAGE:\npython %s folder sample_list name\n"  %os.path.abspath(argv[0]))

############
def multiqc_help():
	"""
	multiQC_ software description message.
	"""
	## [TODO]
	print (colored("\n\n***** TODO: Generate this help message *****\n\n", 'red'))

############
def main():
  	
  	## control if options provided or help
	if len(sys.argv) > 1:
		print ("")
	else:
		help_options()
		exit()    	

	## ARGV
	folder_name = os.path.abspath(argv[1])
	samples_path = os.path.abspath(argv[2])
	name = argv[3]

	## call
	multiQC_call(samples_path, name, folder_name,"")
		
############
'''******************************************'''
if __name__== "__main__":
	main()
	
##Options multiqc, version 1.7 :
##  -f, --force                     Overwrite any existing reports
##  -d, --dirs                      Prepend directory to sample names
##  -dd, --dirs-depth INTEGER       Prepend [INT] directories to sample names. Negative number to take from start of path.
##  -s, --fullnames                 Do not clean the sample names (leave as full file name)
##  -i, --title TEXT                Report title. Printed as page header, used for filename if not otherwise specified.
##  -b, --comment TEXT              Custom comment, will be printed at the top of the report.
##  -n, --filename TEXT             Report filename. Use 'stdout' to print to standard out.
##  -o, --outdir TEXT               Create report in the specified output directory.
##  -t, --template [simple|default_dev|default|sections|geo] Report template to use.
##  --tag TEXT                      Use only modules which tagged with this keyword, eg. RNA
##  --view-tags, --view_tags        View the available tags and which modules they load
##  -x, --ignore TEXT               Ignore analysis files (glob expression)
##  --ignore-samples TEXT           Ignore sample names (glob expression)
##  --ignore-symlinks               Ignore symlinked directories and files
##  --sample-names PATH             File containing alternative sample names
##  -l, --file-list                 Supply a file containing a list of file paths to be searched, one per row
##  -e, --exclude [module name]     Do not use this module. Can specify multiple times.
##  -m, --module [module name]      Use only this module. Can specify multiple times.
##  --data-dir                      Force the parsed data directory to be created.
##  --no-data-dir                   Prevent the parsed data directory from being created.
##  -k, --data-format [json|yaml|tsv] Output parsed data in a different format. Default: tsv
##  -z, --zip-data-dir              Compress the data directory.
##  -p, --export                    Export plots as static images in addition to the report
##  -fp, --flat                     Use only flat plots (static images)
##  -ip, --interactive              Use only interactive plots (HighCharts Javascript)
##  --lint                          Use strict linting (validation) to help code development
##  --pdf                           Creates PDF report with 'simple' template. Requires Pandoc to be installed.
##  --no-megaqc-upload              Don't upload generated report to MegaQC, even if MegaQC options are found
##  -c, --config PATH               Specific config file to load, after those in MultiQC dir / home dir / working dir.
##  --cl-config, --cl_config TEXT   Specify MultiQC config YAML on the command line
##  -v, --verbose                   Increase output verbosity.
##  -q, --quiet                     Only show log warnings
##  --version                       Show the version and exit.
##  -h, --help                      Show this message and exit.

