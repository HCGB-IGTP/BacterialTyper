#usr/bin/env python
'''
This module calls MultiQC to generate HTML statistics reports 
Jose F. Sanchez
Copyright (C) 2019 Lauro Sumoy Lab, IGTP, Spain
'''
## useful imports
import os
import io
import sys
from io import open
from sys import argv
from termcolor import colored

## import my modules
from BacterialTyper import functions
from BacterialTyper import config

############
def multiQC_module_call(givenList, name, path, option):
	pathFile = list2file(givenList, path)
	multiQC_call(pathFile, name, path, option)	
	
############
def multiQC_call(pathFile, name, folder, option):
	multiqc_bin = "multiqc" ## if we activate the environment it should be in $PATH
	## set options for call
	cmd = "%s -o %s -n %s -l %s -p -i 'MultiQC report' -b 'HTML report generated for multiple samples and steps' %s" %(multiqc_bin, folder, name, pathFile, option)
	return(functions.system_call(cmd))

############
def list2file(givenList, path):
	# generate a txt file containing information for 
	name = path + '/' + 'samples.txt'
	outfile = open(name, "w")
	outfile.write("\n".join(givenList))
	outfile.close()
	return(name)

############
def	help_options():
	print ("\nUSAGE:\npython %s folder sample_list name\n"  %os.path.abspath(argv[0]))

############
def multiqc_help():
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
	multiQC_call(samples_path, name, folder_name)
		
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

