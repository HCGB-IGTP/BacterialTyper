#!/usr/bin/env python3
##########################################################
## Jose F. Sanchez										##
## Copyright (C) 2019 Lauro Sumoy Lab, IGTP, Spain		##
##########################################################

## useful imports
import time
import io
import os
import re
import sys
import concurrent.futures
from termcolor import colored
import pandas as pd
import shutil

## import my modules
from BacterialTyper import functions
from BacterialTyper import config
from BacterialTyper import spades_assembler
from BacterialTyper import annotation
from BacterialTyper import BUSCO_caller
from BacterialTyper.modules import qc
from BacterialTyper.modules import help_info

##
global assembly_stats
assembly_stats = {}

####################################
def run_assembly(options):
	"""Main function of the assemble module.
	
	It assembles each sample using SPADES_ and checks quality using BUSCO_ software and database.

	
	.. seealso:: This function depends on other BacterialTyper functions called:
	
		- :func:`BacterialTyper.sampleParser.help_format`
	
		- :func:`BacterialTyper.BUSCO_caller.print_help_BUSCO`
	
		- :func:`BacterialTyper.info.project_help`
		
		- :func:`BacterialTyper.multiQC_report.multiqc_help`
		
		- :func:`BacterialTyper.functions.pipeline_header`
		
		- :func:`BacterialTyper.functions.boxymcboxface`
		
		- :func:`BacterialTyper.functions.print_time`
		
		- :func:`BacterialTyper.sampleParser.get_files`
		
		- :func:`BacterialTyper.functions.outdir_project`
		
		- :func:`BacterialTyper.functions.optimize_threads`
		
		- :func:`BacterialTyper.modules.assemble.check_sample_assembly`
		
		- :func:`BacterialTyper.functions.create_subfolder`
		
		- :func:`BacterialTyper.modules.qc.BUSCO_check`
		
	.. include:: ../../links.inc	 	
	
	"""
	
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
		help_info.project_help()
		exit()
	elif (options.help_multiqc):
		## information for Multiqc
		multiQC_report.multiqc_help()
		exit()
	
	### set as default paired_end mode
	if (options.single_end):
		options.pair = False
	else:
		options.pair = True

	## message header
	functions.pipeline_header()
	functions.boxymcboxface("Assembly module")
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

	## get files
	pd_samples_retrieved = sampleParser.get_files(options, input_dir, "trim", ['_trim_'])
	
	## debug message
	if (Debug):
		print (colored("**DEBUG: pd_samples_retrieve **", 'yellow'))
		print (pd_samples_retrieved)
	
	## generate output folder, if necessary
	print ("\n+ Create output folder(s):")
	if not options.project:
		functions.create_folder(outdir)
	outdir_dict = functions.outdir_project(outdir, options.project, pd_samples_retrieved, "assemble")

	### call assemble using spades
	start_time_partial = start_time_total
	start_time_partial_assembly = start_time_partial
	
	## optimize threads
	name_list = set(pd_samples_retrieved["name"].tolist())
	threads_job = functions.optimize_threads(options.threads, len(name_list)) ## threads optimization
	max_workers_int = int(options.threads/threads_job)

	## debug message
	if (Debug):
		print (colored("**DEBUG: options.threads " +  str(options.threads) + " **", 'yellow'))
		print (colored("**DEBUG: max_workers " +  str(max_workers_int) + " **", 'yellow'))
		print (colored("**DEBUG: cpu_here " +  str(threads_job) + " **", 'yellow'))

	# Group dataframe by sample name
	sample_frame = pd_samples_retrieved.groupby(["name"])

	# We can use a with statement to ensure threads are cleaned up promptly
	print ('+ Running modules SPADES...')
	with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers_int) as executor:
		## send for each sample
		commandsSent = { executor.submit( check_sample_assembly, name, outdir_dict[name],  sorted(cluster["sample"].tolist()), threads_job): name for name, cluster in sample_frame }

		for cmd2 in concurrent.futures.as_completed(commandsSent):
			details = commandsSent[cmd2]
			try:
				data = cmd2.result()
			except Exception as exc:
				print ('***ERROR:')
				print (cmd2)
				print('%r generated an exception: %s' % (details, exc))
		
	## functions.timestamp
	print ("\n+ Assembly of all samples finished: ")
	start_time_partial = functions.timestamp(start_time_partial_assembly)
	
	##
	if (assembly_stats):
		## get all assembly stats
		outdir_report = functions.create_subfolder("report", outdir)
		final_dir = functions.create_subfolder("assembly_stats", outdir_report)
		final_sub_dir = functions.create_subfolder("samples", final_dir)
		

		#### summary and information
		results_summary_toPrint_all = pd.DataFrame()		
		column_names = ("Total Sequences", "Total Length (bp)", "Total length (no Ns)", 
                "Adenine (A)", "Thymine (T)", "Cytosine (C)", "Guanine (G)", "A+T", "C+G", "Any Nucleotide (N)", 
                "Capture Gaps", "Capture Gaps Length", "Capture Gaps Length/Total Length (%)", 
                "Number Seqs > 10 kb", "% Seqs > 10 kb", "Total Length (bp) > 10 kb", "% Bases > 10 kb", 
                "MinLen", "MaxLen", "Average Len", "Median Len", "N50", "L50", "Length (no Ns)")
		
		for stats in assembly_stats:
			##
			file_stats = assembly_stats[stats]
			shutil.copy(file_stats, final_sub_dir)
		
			tmp_name = os.path.splitext(file_stats)
			#csv
			csv_file = tmp_name[0] + '.csv'
		
			## subset dataframe	& print result
			results_summary_toPrint = pd.read_csv(csv_file, sep=',', header=None).transpose()
			results_summary_toPrint.columns = column_names
			results_summary_toPrint['sample'] = stats
		
			## add all data
			results_summary_toPrint_all = pd.concat([results_summary_toPrint_all, results_summary_toPrint], ignore_index=True)

		## write to excel
		name_excel_summary = final_dir + '/summary_stats.xlsx'
		writer_summary = pd.ExcelWriter(name_excel_summary, engine='xlsxwriter') ## open excel handle
		
		## filter important columns		
		results_summary_toPrint_all = results_summary_toPrint_all.set_index('sample')			
		results_summary_toPrint_all = results_summary_toPrint_all[["Total Sequences", "Total Length (bp)", "Total length (no Ns)",
														"Adenine (A)", "Thymine (T)", "Cytosine (C)", "Guanine (G)", "A+T", "C+G", "Any Nucleotide (N)",
                                                  		"MinLen", "MaxLen", "Average Len", "Median Len", "N50", "L50",
                                                    	"Number Seqs > 10 kb", "% Seqs > 10 kb", "Total Length (bp) > 10 kb", "% Bases > 10 kb"]]
		
		
		## >10 kb subset
		subset_summary_toPrint_all = results_summary_toPrint_all[["Total Sequences", "Total Length (bp)", "Number Seqs > 10 kb", "% Bases > 10 kb", "Median Len", "N50", "L50"]]
		subset_summary_toPrint_all.to_excel(writer_summary, sheet_name="summary") ## write excel handle
		
		## nucleotide tab
		nucleotides_summary_toPrint_all = results_summary_toPrint_all[["Adenine (A)", "Thymine (T)", "Cytosine (C)", "Guanine (G)", "A+T", "C+G", "Any Nucleotide (N)"]]
		nucleotides_summary_toPrint_all.to_excel(writer_summary, sheet_name="nucleotides") ## write excel handle
		
		results_summary_toPrint_all.to_excel(writer_summary, sheet_name="all_data") ## write excel handle
		
		writer_summary.save() ## close excel handle
	
	### symbolic links
	print ("+ Retrieve all genomes assembled...")
	
	### BUSCO check assembly
	results = qc.BUSCO_check(outdir, outdir, options, start_time_partial, "genome")
	
	## print to file results	 	
	print ("\n*************** Finish *******************")
	start_time_partial = functions.timestamp(start_time_total)

	print ("+ Exiting Assembly module.")
	exit()

#############################################
def check_sample_assembly(name, sample_folder, files, threads):
	"""Checks if sample is assembled.
	
	It checks wether a sample is assembled or not by reading file *sample_folder/.success_all* (generated by :func:`BacterialTyper.functions.print_time_stamp`). 
	
	If file not available (no previous assembly or not suceeded it) it calls :func:`BacterialTyper.spades_assembler.run_module_assembly` to generate assembly for the sample speficied.
	
	:param name: Sample name or tag to identify sample.
	:param sample_folder:  directory to generate assembly ouptut. It must exist.
	:param files: List containing files (fastq R1 & R2) for the sample to be assembled.
	:param threads: Number of CPUs to use
	:type name: string
	:type sample_folder: string 
	:type files: list
	:type threads: integer
	
	:return: Populates dictionary assembly_stats with assembly stats information file
	:rtype: Dataframe
	
	.. seealso:: This function depends on other BacterialTyper functions called:
	
		- :func:`BacterialTyper.spades_assembler.run_module_assembly`
	
		- :func:`BacterialTyper.functions.print_time_stamp`
	
	"""
	## check if previously assembled and succeeded
	filename_stamp = sample_folder + '/.success_all'
	if os.path.isfile(filename_stamp):
		stamp =	functions.read_time_stamp(filename_stamp)
		print (colored("\tA previous command generated results on: %s [%s]" %(stamp, name), 'yellow'))
		assembly_stats[name] = sample_folder + '/' + name + "_assembly.fna_stats.txt"

	else:
	
		## debug message
		if (Debug):
			print (colored("**DEBUG: spades_assembler.run_module_assembly call**", 'yellow'))
			print ("spades_assembler.run_module_assembly(name, sample_folder, files[0], files[1], threads)")
			print ("spades_assembler.run_module_assembly " + name + "\t" + sample_folder + "\t" + files[0] + "\t" + files[1] + "\t" +str(threads) + "\n")
	
		# Call spades_assembler
		code = spades_assembler.run_module_assembly(name, sample_folder, files[0], files[1], threads)
		
		if (code != 'FAIL'):
			## success stamps
			filename_stamp = sample_folder + '/.success_all'
			stamp =	functions.print_time_stamp(filename_stamp)
			assembly_stats[name] = code

