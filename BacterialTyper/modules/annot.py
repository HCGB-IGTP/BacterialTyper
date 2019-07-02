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

	### species_identification_KMA -> most similar taxa
	functions.pipeline_header()
	functions.boxymcboxface("Assembly annotation")

	print ("--------- Starting Process ---------")
	functions.print_time()

	## print further information for Prokka	
	if (options.help_Prokka):
		annotation.print_list_prokka()
		exit()
	
	## print further information for BUSCO databases	
	if (options.help_BUSCO):
		BUSCO_caller.print_help_BUSCO()
		exit()
	
	## absolute path for in & out
	input_dir = os.path.abspath(options.input)
	outdir = os.path.abspath(options.output_folder)

	## generate output folder
	functions.create_folder(outdir)

	### time stamp
	start_time_partial = start_time_total
	
	### symbolic links
	print ("+ Retrieve all genomes assembled...")
	
	### batch provided
	if options.batch:
		## csv file containing sample name and file path
		pd_samples_retrieved = pd.read_csv(options.batch, sep=',',header=None)
		pd_samples_retrieved.columns = ["samples", "tag", "file"]
	else:
		fasta_ext = ('.fna', '.fasta')
		genome_pd_samples_retrieved = sample_prepare.get_files(options, input_dir, "chromosome", fasta_ext)
	
		print ("+ Retrieve all plasmids assembled...")
		plasmid_pd_samples_retrieved = sample_prepare.get_files(options, input_dir, "plasmid", fasta_ext)
		
		frames = [plasmid_pd_samples_retrieved, genome_pd_samples_retrieved]
		pd_samples_retrieved = pd.concat(frames, ignore_index=True)

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
	threads_module = functions.optimize_threads(options.threads, pd_samples_retrieved.index.size)

	## send for each sample
	with concurrent.futures.ThreadPoolExecutor(max_workers=int(options.threads)) as executor:
		commandsSent = { executor.submit(annotation.module_call, row['file'], options.kingdom, options.genera, outdir + '/' + row['samples'] + '_' + row['tag'], row['samples'] , threads_module): index for index, row in pd_samples_retrieved.iterrows() }
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

	## retrieve information
	givenList = []
	protein_files = []
	print ("+ Detail information for each sample could be identified in separate folders:")
	for index, row in pd_samples_retrieved.iterrows():
		fold_sample = outdir + '/' + row['samples'] + '_' + row['tag'] + '/'
		givenList.append(fold_sample)
		if row['tag'] != 'plasmid':
			protein_files.extend(functions.retrieve_matching_files(fold_sample, '.faa'))
		print ('\t- %s' %fold_sample)	

	print (protein_files)

	### report generation
	if (options.skip_report):
		print ("+ No annotation report generation...")
	else:
		### report generation
		functions.boxymcboxface("Annotation report")

		outdir_report = functions.create_subfolder("report", outdir)
		
		## check if previously report generated
		filename_stamp = outdir_report + '/.success'
		done=0
		if os.path.isdir(outdir_report):
			if os.path.isfile(filename_stamp):
				stamp =	functions.read_time_stamp(filename_stamp)
				print (colored("\tA previous report generated results on: %s" %stamp, 'yellow'))
				done=1
		
		## generate report
		if done==0:
			## get subdirs generated and call multiQC report module
			multiQC_report.multiQC_module_call(givenList, "Prokka", outdir_report, "-dd 1")
			print ('\n+ A summary HTML report of each sample is generated in folder: %s' %outdir_report)
		
			## success stamps
			filename_stamp = outdir_report + '/.success'
			stamp =	functions.print_time_stamp(filename_stamp)

	## time stamp
	start_time_partial_BUSCO = functions.timestamp(start_time_total)

	## Check each annotation using BUSCO
	functions.boxymcboxface("BUSCO Annotation Quality check")
	database_folder = os.path.abspath(options.database)
	outdir_BUSCO = functions.create_subfolder("BUSCO", outdir)
	BUSCO_Database = database_folder + '/BUSCO'
	dataFrame_results = qc.BUSCO_call(options.BUSCO_dbs, protein_files, BUSCO_Database, outdir_BUSCO, options.threads, "proteins")
	
	## functions.timestamp
	print ("+ Quality control of all samples finished: ")
	start_time_partial = functions.timestamp(start_time_partial_BUSCO)
	
	## summarize
	if (options.skip_report):
		print ("+ No report generation...")
	else:
		## generate plots
		print ("+ Generate summarizing plots...")
		qc.BUSCO_plots(dataFrame_results, outdir, options.threads)	
	
	## multiqc report plot
	## busco plot of samples or datasets per sample

	print ("\n*************** Finish *******************")
	start_time_partial = functions.timestamp(start_time_total)

	print ("+ Exiting Annotation module.")
	exit()

