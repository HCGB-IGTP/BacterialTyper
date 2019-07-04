#!/usr/bin/env python3
'''
This code calls Prokka to generate a functional assembly annotation
Jose F. Sanchez
Copyright (C) 2019 Lauro Sumoy Lab, IGTP, Spain
'''
## useful imports
import time
import io
import os
import re
import sys
import concurrent.futures
from termcolor import colored
import pandas as pd

## import my modules
from BacterialTyper import functions
from BacterialTyper import config
from BacterialTyper import annotation
from BacterialTyper import BUSCO_caller
from BacterialTyper.modules import qc
from BacterialTyper import multiQC_report
from BacterialTyper.modules import sample_prepare

####################################
def run(options):

	## init time
	start_time_total = time.time()

	## debugging messages
	global Debug
	if (options.debug):
		Debug = True
	else:
		Debug = False

	##################################
	### show help messages if desired	
	##################################
	if (options.help_format):
		## help_format option
		sampleParser.help_format()
		exit()

	elif (options.help_BUSCO):
		## information for BUSCO
		BUSCO_caller.print_help_BUSCO()
		exit()

	elif (options.help_project):
		## information for project
		help.project_help()
		exit()

	elif (options.help_multiqc):
		## information for Multiqc
		multiQC_report.multiqc_help()

	elif (options.help_Prokka):
		## information for Prokka
		annotation.print_list_prokka()
		exit()

	### 
	functions.pipeline_header()
	functions.boxymcboxface("Assembly annotation")

	print ("--------- Starting Process ---------")
	functions.print_time()

	## absolute path for in & out
	input_dir = os.path.abspath(options.input)
	outdir=""

	## set mode: project/detached
	if (options.project):
		outdir = input_dir		
	elif (options.detached):
		outdir = os.path.abspath(options.output_folder)

	### symbolic links
	print ("+ Retrieve all genomes assembled...")

	## get files
	pd_samples_retrieved = sample_prepare.get_files(options, input_dir, "assembly", "fna")

	## debug message
	if (Debug):
		print (colored("**DEBUG: pd_samples_retrieve **", 'yellow'))
		print (pd_samples_retrieved)
	
	## generate output folder, if necessary
	print ("\n+ Create output folder(s):")
	outdir_dict = functions.outdir_project(outdir, options.project, pd_samples_retrieved, "annot")
	
	if not options.project:
		functions.create_folder(outdir)

	## annotate
	print ("+ Annotate assemblies using prokka:")
	print ("\t-Option: kingdom = ", options.kingdom,"; Annotation mode")
	if options.genera == 'Other':
		print ("\t-Option: genera = Off; No genus-specific BLAST databases option provided")
	else:
		print ("\t-Option: genera = ", options.genera,"; Genus-specific BLAST databases option provided")

	print ("\t-Option: addgenes; Add 'gene' features for each 'CDS' feature")
	print ("\t-Option: addmrna;  Add 'mRNA' features for each 'CDS' feature")
	print ("\t-Option: cdsrnaolap;  Allow [tr]RNA to overlap CDS")

	## optimize threads
	## workers:
	name_list = set(pd_samples_retrieved["name"].tolist())
	max_workers_int = len(name_list)
	
	## number_samples = pd_samples_retrieved.index.size => Number samples
	threads_module = functions.optimize_threads(options.threads, max_workers_int) ## fix threads_module optimization

	## debug message
	if (Debug):
		print (colored("**DEBUG: options.threads " +  str(options.threads) + " **", 'yellow'))
		print (colored("**DEBUG: max workers " +  str(max_workers_int) + " **", 'yellow'))
		print (colored("**DEBUG: cpu/process " +  str(threads_module ) + " **", 'yellow'))

	## send for each sample
	with concurrent.futures.ThreadPoolExecutor(max_workers=int(max_workers_int)) as executor:
		commandsSent = { executor.submit(annot_caller, row['sample'], outdir_dict[row['name']], options, row['name'], threads_module): index for index, row in pd_samples_retrieved.iterrows() }
		for cmd2 in concurrent.futures.as_completed(commandsSent):
			details = commandsSent[cmd2]
			try:
				data = cmd2.result()
			except Exception as exc:
				print ('***ERROR:')
				print (cmd2)
				print('%r generated an exception: %s' % (details, exc))

	## time stamp
	start_time_partial = functions.timestamp(start_time_total)

	## get folders
	givenList = [ v for v in outdir_dict.values() ]
	protein_files = []
	print ("+ Detail information for each sample could be identified in separate folders:")
	for folder in givenList:
		print ('\t + ', folder)
		protein_files.extend(functions.retrieve_matching_files(folder, '.faa'))

	### report generation
	if (options.skip_report):
		print ("+ No annotation report generation...")
	else:
		### report generation
		functions.boxymcboxface("Annotation report")
		outdir_report = functions.create_subfolder("report", outdir)
		
		PROKKA_report = functions.create_subfolder("annotation", outdir_report)
		print ('\n+ A summary HTML report of each sample is generated in folder: %s' %PROKKA_report)
		
		## check if previously report generated
		filename_stamp = PROKKA_report + '/.success'
		done=0
		if os.path.isdir(PROKKA_report):
			if os.path.isfile(filename_stamp):
				stamp =	functions.read_time_stamp(filename_stamp)
				print (colored("\tA previous report generated results on: %s" %stamp, 'yellow'))
				done=1
		
		## generate report
		if done==0:
			## get subdirs generated and call multiQC report module
			multiQC_report.multiQC_module_call(givenList, "Prokka", PROKKA_report, "-dd 2")
			print ('\n+ A summary HTML report of each sample is generated in folder: %s' %PROKKA_report)
		
			## success stamps
			filename_stamp = PROKKA_report + '/.success'
			stamp =	functions.print_time_stamp(filename_stamp)

	## time stamp
	start_time_partial_BUSCO = functions.timestamp(start_time_total)

	## Check each annotation using BUSCO
	results = qc.BUSCO_check(input_dir, outdir, options, start_time_partial_BUSCO, "proteins")

	## print to file: results	 
	
	print ("\n*************** Finish *******************")
	start_time_partial = functions.timestamp(start_time_total)

	print ("+ Exiting Annotation module.")
	exit()


#############################################
def annot_caller(seq_file, sample_folder, options, name, threads):
	## check if previously assembled and succeeded
	filename_stamp = sample_folder + '/.success'

	if os.path.isfile(filename_stamp):
		stamp =	functions.read_time_stamp(filename_stamp)
		print (colored("\tA previous command generated results on: %s" %stamp, 'yellow'))
	else:
	
		## debug message
		if (Debug):
			print (colored("**DEBUG: annotation.module_call call**", 'yellow'))
			print (" annotation.module_call (seq_file, options.kingdom, options.genera, sample_folder, name, threads)")
			print (" annotation.module_call " + seq_file + "\t" + options.kingdom + "\t" + options.genera + "\t" + sample_folder + "\t" + name + "\t" + str(threads))

		# Call annotation
		annotation.module_call(seq_file, options.kingdom, options.genera, sample_folder, name, threads)
		
		


