#!/usr/bin/env python3
##############################################################
## Jose F. Sanchez											##
## Copyright (C) 2019-2020 Lauro Sumoy Lab, IGTP, Spain		##
##############################################################
"""
This code calls several programs to generate a Mobile Genetic Element analysis
"""

## useful imports
import time
import io
import os
import re
import sys
from io import open
import concurrent.futures
from termcolor import colored
import pandas as pd

## import my modules
from BacterialTyper import functions, bacteriophage
from BacterialTyper import config
from BacterialTyper import annotation
from BacterialTyper import bacteriophage
from BacterialTyper import sampleParser
from BacterialTyper.modules import help_info
from BacterialTyper import database_generator
from BacterialTyper import database_user

####################################
def run_MGE(options):

	##################################
	### show help messages if desired	
	##################################
	if (options.help_format):
		## help_format option
		sampleParser.help_format()
		exit()

	if (options.help_project):
		## information for project
		help_info.project_help()
		exit()
	
	if (options.help_PhiSpy):
		## information for PhiSpy software
		bacteriophage.help_PhiSpy()
		exit()

	if (options.help_MGE_analysis):
		## information for MGE module analysis
		help_MGE_analysis()
		exit()

	if (options.help_input_MGE):
		## information for PhiSpy
		help_input_MGE()
		exit()

	## debugging messages
	global Debug
	if (options.debug):
		Debug = True
	else:
		Debug = False
		
	### set as default paired_end mode
	if (options.single_end):
		options.pair = False
	else:
		options.pair = True

	## init time
	start_time_total = time.time()

	### species_identification_KMA -> most similar taxa
	functions.pipeline_header()
	functions.boxymcboxface("Mobile Genetic Elements (MGE) identification")

	print ("--------- Starting Process ---------")
	functions.print_time()
	
	## absolute path for in & out
	input_dir = os.path.abspath(options.input)
	outdir=""

	### init dataframe
	pd_samples_retrieved = pd.DataFrame()

	## set mode: project/detached
	global Project
	if (options.project):
		outdir = input_dir		
		Project=True

		## get user data	
		pd_samples_retrieved = database_user.get_userData_files(options, input_dir)
		## returns merge dafarame 

	elif (options.detached):
		#########
		Project=False
		outdir = os.path.abspath(options.output_folder)
		## input bactch
		if (os.path.isfile(input_dir)):			
			input_data = pd.read_csv(input_dir, header=None, sep=',')
			input_data.columns = ["name", "ext", "sample"]

			## debug message
			if (Debug):
				print (colored("**DEBUG: input_data csv **", 'yellow'))
				print (input_data)
			
			## get info from reads			
			reads_dataframe = input_data[input_data['ext']=='reads']
			if not reads_dataframe.empty:
				file_name_list = reads_dataframe['sample'].tolist()
				pd_samples_reads = sampleParser.get_fields(file_name_list, options.pair, Debug)

				## discard lane information
				pd_samples_reads = pd_samples_reads.drop(["lane", "lane_file", "gz"], axis=1)

				## debug message
				if (Debug):
					print (colored("**DEBUG: reads_dataframe **", 'yellow'))
					print (reads_dataframe)
					print (colored("**DEBUG: pd_samples_reads **", 'yellow'))
					print (pd_samples_reads)
			## init to 0
			else: 
				pd_samples_reads = pd.DataFrame()

			## other data
			other_dataframe = input_data[input_data['ext']!='reads']
			if not other_dataframe.empty:
			
				## assembly data
				pd_samples_assembly = other_dataframe[other_dataframe['ext'] == 'assembly'] ## [TODO: Fix SettingWithCopyWarning]

				if not pd_samples_assembly.empty:
					pd_samples_assembly = pd_samples_assembly.rename(index=str, columns={'ext':'tag'}) 
					pd_samples_assembly.loc[:,'ext'] = 'fna'
					pd_samples_assembly.loc[:,'dirname'] = pd_samples_assembly.apply(lambda row: os.path.dirname(row['sample']), axis=1)
			
					## debug message
					if (Debug):
						print (colored("**DEBUG: pd_samples_assembly **", 'yellow'))
						print (pd_samples_assembly)

				## init to 0
				else: 
					pd_samples_assembly = pd.DataFrame()

				## annot_data
				pd_samples_annot = other_dataframe[other_dataframe['ext']!='assembly']  ## [TODO: Fix SettingWithCopyWarning]
				
				if not pd_samples_annot.empty:
					pd_samples_annot.loc[:,'tag'] = 'annot'
					pd_samples_annot.loc[:,'dirname'] = pd_samples_annot.apply(lambda row: os.path.dirname(row['sample']), axis=1)
					
					## debug message
					if (Debug):
						print (colored("**DEBUG: pd_samples_annot **", 'yellow'))
						print (pd_samples_annot)

				## init to 0
				else: 
					pd_samples_annot = pd.DataFrame()

			## init to 0
			else: 
				pd_samples_assembly = pd.DataFrame()
				pd_samples_annot = pd.DataFrame()

		else:
			print (colored("\n\n**** ERROR: input provided is not a readable file. Please provide a csv if --detached mode. See --help_input_MGE for further details", 'red'))
			exit()


		## merge into pd_samples_retrieved
		frames = [pd_samples_reads, pd_samples_assembly, pd_samples_annot]
		pd_samples_retrieved = pd.concat(frames, sort=True, join='outer')
		
	########

	## debug message
	if (Debug):
		print (colored("**DEBUG: pd_samples_retrieved **", 'yellow'))
		print (pd_samples_retrieved)
		
	## generate output folder, if necessary
	print ("\n+ Create output folder(s):")
	if not options.project:
		functions.create_folder(outdir)

	## for each sample
	outdir_dict = functions.outdir_project(outdir, options.project, pd_samples_retrieved, "MGE")
	
	### time stamp
	start_time_partial = start_time_total
	start_time_partial_assembly = start_time_partial
	
	## once all data is retrieved... 
	print ("\n+ Mobile genetic elements (MGE) analysis:")	
	global plasmid_bool
	global bacteriophage_bool
	global genomic_island_bool

	## perform different analysis
	if (options.plasmid_analysis):
		print ("\t- Option: Plasmid analysis")
		plasmid_bool = True
		
	elif (options.phage_analysis):
		bacteriophage_bool = True
		print ("\t- Option: Phage analysis")
		
	elif (options.GI_analysis):
		genomic_island_bool = True
		print ("\t- Option: Genomic Island analysis")

	elif (options.all_data):
		print ("\t- Option: Plasmid analysis")
		print ("\t- Option: Phage analysis")
		print ("\t- Option: Genomic Island analysis")
		
		plasmid_bool = True
		bacteriophage_bool = True
		genomic_island_bool = True

	## parameters 
	print ("+ Parameters:")
	if (bacteriophage_bool):
		print ("\tPhiSpy Training set: ", str(options.training_set))
		print ("\tPhiSpy Window size: ", str(options.window_size))
		print ("\tPhiSpy Phage genes: ", str(options.phage_genes))
		
	print ("")
	print ("")

	## re-index dataframe
	pd_samples_retrieved = pd_samples_retrieved.reset_index()
	
	## debug message
	if (Debug):
		print (colored("**DEBUG: pd_samples_retrieved **", 'yellow'))
		print (pd_samples_retrieved)
			
	## optimize threads
	name_list = set(pd_samples_retrieved["name"].tolist())
	threads_job = functions.optimize_threads(options.threads, len(name_list)) ## threads optimization
	max_workers_int = int(options.threads/threads_job)

	# Group dataframe by sample name	
	sample_frame = pd_samples_retrieved.groupby(["name"])

	## debug message
	if (Debug):
		print (colored("**DEBUG: options.threads " +  str(options.threads) + " **", 'yellow'))
		print (colored("**DEBUG: max_workers " +  str(max_workers_int) + " **", 'yellow'))
		print (colored("**DEBUG: cpu_here " +  str(threads_job) + " **", 'yellow'))

	## send for each sample
	with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers_int) as executor:
		commandsSent = { executor.submit(MGE_caller, outdir_dict[name], name, options, threads_job, cluster): name for name, cluster in sample_frame }
		for cmd2 in concurrent.futures.as_completed(commandsSent):
			details = commandsSent[cmd2]
			try:
				data = cmd2.result()
			except Exception as exc:
				print ('***ERROR:')
				print (cmd2)
				print('%r generated an exception: %s' % (details, exc))

		exit()

###########################
def MGE_caller(output_dir, name, options, threads, dataFrame_sample):
	"""Identify Mobile Genetic Elements
	
	This function is called for each sample to assess the putative plasmids, inserted phages and genomic islands.
	
	.. seealso: This functions depend on third party software called using several BacterialTyper functions such as:
	
		- :func:`BacterialTyper.bacteriophage.call_ident_bacteriophage`
		
		- :func:`BacterialTyper.bacteriophage.results_PhiSpy`
		
		- :func:`BacterialTyper.`
		
		- :func:`BacterialTyper.`
	
	"""
	
	print ("+ MGE analysis for sample: ", name)
	print ("")
	
	if plasmid_bool:

		functions.print_sepLine("*",50, False)
		print (' Plasmid identification analysis')
		functions.print_sepLine("*",50, False)

		## for each sample
		outdir_dict_plasmid = functions.outdir_subproject(output_dir, dataFrame_sample, "plasmid")

		## debug message
		if (Debug):
			print (colored("**DEBUG: Dir" + str(outdir_dict_plasmid), 'yellow'))
			
			## ToDo implement plasmid analysis
			
			
		print ("")
		
	if bacteriophage_bool:
		
		functions.print_sepLine("*",50, False)
		print (' Phage identification analysis')
		functions.print_sepLine("*",50, False)
		
		## for each sample
		outdir_dict_phage = functions.outdir_subproject(output_dir, dataFrame_sample, "phage")

		## debug message
		if (Debug):
			print (colored("**DEBUG: Dir" + str(outdir_dict_phage), 'yellow'))

		## get Genbank file generated with PROKKA
		gbk_file = dataFrame_sample.loc[dataFrame_sample['ext'] == 'gbf']['sample'].tolist() ## [TODO: Fix SettingWithCopyWarning]
		
		## Call phispy
		bacteriophage.ident_bacteriophage(gbk_file[0], name, outdir_dict_phage[name], options.training_set, Debug, 
											window_size=options.window_size, number_phage_genes=options.phage_genes)
		
		## Parse results
		bacteriophage.results_PhiSpy(outdir_dict_phage[name], name)

	if genomic_island_bool:

		functions.print_sepLine("*",50, False)
		print (' Genomic island identification analysis')
		functions.print_sepLine("*",50, False)

		## for each sample
		outdir_dict_GI = functions.outdir_subproject(output_dir, dataFrame_sample, "genomic_island")

		## debug message
		if (Debug):
			print (colored("**DEBUG: Dir"+ str(outdir_dict_GI), 'yellow'))
		
		
			## ToDo implement genomic_island analysis
		
		print ("")

###########################
def help_MGE_analysis():
	print (colored("\n\n***** TODO: Generate this help message *****\n\n", 'red'))

###########################
def help_input_MGE():
	message = '''
##################
  Project mode
##################
Project a folder containing the project previously generated using the different modules: prep, trimm, assemble, annot, ident, profile...

##################
  Detached mode: 
##################
Provide a csv file containing the sample name, the type of file and the path for each file as input. Set --batch option. 
Provide as type of file: reads, gff, gbf, faa, assembly. 
Do not provide any header.  

e.g.: 

---
sample1,gff,/path/to/sample1/annotation/file.gff
sample1,reads,/path/to/sample1/reads/name_R1.fastq
sample1,reads,/path/to/sample1/reads/name_R2.fastq
sample1,assembly,/path/sample1/name.fasta
sample1,faa,/path/sample1/annotation/name.faa
sample1,gbf,/path/sample1/annotation/name.gbf
sample2,gff,/path/to/sample2/annotation/file.gff
sample2,reads,/path/to/sample2/reads/name_R1.fastq
sample2,reads,/path/to/sample2/reads/name_R2.fastq
sample2,assembly,/path/sample2/name.fasta
sample2,faa,/path/sample2/annotation/name.faa
sample2,gbf,/path/sample2/annotation/name.gbf
---
'''

	print (message)

