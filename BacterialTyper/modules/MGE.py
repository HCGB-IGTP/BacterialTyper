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
import shutil

## import my modules
from BacterialTyper.scripts import bacteriophage, genomic_island
from BacterialTyper.config import set_config
from BacterialTyper.scripts import annotation
from BacterialTyper.scripts import bacteriophage
from BacterialTyper.modules import help_info
from BacterialTyper.scripts import database_generator
from BacterialTyper.scripts import database_user
from BacterialTyper import __version__ as pipeline_version

###
global phage_Results
phage_Results = {}
global GI_Results
GI_Results = {}


####################################
def run_MGE(options):

	##################################
	### show help messages if desired	
	##################################
	if (options.help_format):
		## help_format option
		help_info.help_fastq_format()
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

	if (options.help_Dimob):
		## information for PhiSpy
		help_input_Dimob()
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
	HCGB_aes.pipeline_header("BacterialTyper", ver=pipeline_version)
	HCGB_aes.boxymcboxface("Mobile Genetic Elements (MGE) identification")

	print ("--------- Starting Process ---------")
	functions.print_time()
	
	## absolute path for in & out
	input_dir = os.path.abspath(options.input)
	outdir=""

	### init dataframe
	pd_samples_retrieved = pd.DataFrame()

	## Project mode as default
	if (options.detached):
		options.project = False		
	else:
		options.project = True
		
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
	outdir_dict = functions.outdir_project(outdir, options.project, pd_samples_retrieved, "MGE", options.debug)
	
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
		print ("+ PhiSpy:")
		print ("\tPhiSpy Training set: ", str(options.training_set))
		print ("\tPhiSpy Window size: ", str(options.window_size))
		print ("\tPhiSpy Phage genes: ", str(options.phage_genes))
		
	if (genomic_island_bool):
		print ("+ Dimob:")
		print ("\tDimob dinucleotide bias cutoff (%):", str(options.cutoff_dinuc_bias))
		print ("\tDimob Genomic Island minimum length (bp):", str(options.min_length))
		
	if (plasmid_bool):
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
	max_workers_int = 8

	## there is a problem with RAM, we would set one sample at a time until satisfied

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

	## create summary in report for each
	## get all analysis stats
	print ("+ Create summary of all samples analyzed.")
	outdir_report = functions.create_subfolder("report", outdir)
	final_dir = functions.create_subfolder("MGE_analysis", outdir_report)
	
	if Debug:
		print ("*** DEBUG: MGE_Results ***")
		print (phage_Results)
		print (GI_Results)
	
	results_summary = pd.DataFrame(columns=("sample", "phage", "genomic_island", "plasmid"))
	results_summary = results_summary.set_index("sample")
	
	if (bacteriophage_bool):
		
		phages_dir = functions.create_subfolder("phages", final_dir)
		print ("+ See details of phages results in folder:")
		print ("\t", phages_dir)
		
		for key,value in phage_Results.items():
			file_coord=os.path.join(value, key + '_PhiSpy-prophage-coordinates.tsv')
			if functions.is_non_zero_file(file_coord):
				shutil.copy(file_coord, phages_dir)
				coordinates_pd = pd.read_csv(file_coord, sep="\t")
				results_summary.loc[key, "phage"] = len(coordinates_pd)
	
	if (genomic_island_bool):
		GI_dir = functions.create_subfolder("genomicIslands", final_dir)
		print ("+ See details of genomic islands results in folder:")
		print ("\t", GI_dir)
		
		for key,value in GI_Results.items():
			file_coord=os.path.join(value, key + '.txt')
			if functions.is_non_zero_file(file_coord):
				shutil.copy(file_coord, GI_dir)
				coordinates_pd = pd.read_csv(file_coord, sep="\t")
				results_summary.loc[key, "genomic_island"] = len(coordinates_pd)
	
	if (plasmid_bool):
		plasmid_dir = functions.create_subfolder("plasmids", final_dir)
		print ("+ See details of plasmids results in folder:")
		print ("\t", plasmid_dir)
		
	## use BioCircos to represent this information?
	## Example: 
	## https://github.com/molevol-ub/BacterialDuplicates/blob/master/scripts/R/BioCircos_plotter.R	
	
	print ("+ Printing summary results in XLSX file")
	print ("+ See details in directory: ", final_dir)
	name_excel_summary = os.path.join(final_dir, "MGE_summary.xlsx")
	writer_summary = pd.ExcelWriter(name_excel_summary, engine='xlsxwriter') ## open excel handle
	results_summary.to_excel(writer_summary, sheet_name="summary") ## write excel handle
	writer_summary.save() ## close excel handle

	return()

###########################
def MGE_caller(output_dir, name, options, threads, dataFrame_sample):
	"""Identify Mobile Genetic Elements
	
	This function is called for each sample to assess the putative plasmids, inserted phages and genomic islands.
	
	.. seealso: This functions depend on third party software called using several BacterialTyper functions such as:
	
		- :func:`BacterialTyper.scripts.bacteriophage.call_ident_bacteriophage`
		
		- :func:`BacterialTyper.scripts.bacteriophage.results_PhiSpy`
		
		- :func:`BacterialTyper.scripts.`
		
		- :func:`BacterialTyper.scripts.`
	
	"""
	
	#print ("+ MGE analysis for sample: ", name)
	
	## get Genbank file generated with PROKKA
	gbk_file = dataFrame_sample.loc[dataFrame_sample['ext'] == 'gbf']['sample'].tolist() ## [TODO: Fix SettingWithCopyWarning]
	
	####################
	## plasmid analysis
	####################
	if plasmid_bool:

		#functions.print_sepLine("*",50, False)
		#print ('Plasmid identification analysis')
		#functions.print_sepLine("*",50, False)

		## for each sample
		outdir_plasmid = functions.create_subfolder("plasmid", output_dir)

		## debug message
		if (Debug):
			print (colored("**DEBUG: Dir: " + str(outdir_plasmid), 'yellow'))
			
			## ToDo implement plasmid analysis
			
			
		print ("")
		
	####################
	## phage analysis
	####################
	if bacteriophage_bool:
		
		#functions.print_sepLine("*",50, False)
		#print (' Phage identification analysis')
		#functions.print_sepLine("*",50, False)
		
		## for each sample
		outdir_phage = functions.create_subfolder("phage", output_dir)

		## debug message
		if (Debug):
			print (colored("**DEBUG: Dir: " + str(outdir_phage), 'yellow'))
	
		##
		filename_stamp = outdir_phage + '/.PhiSpy_results'
		# check if previously done
		if os.path.isfile(filename_stamp):
			stamp =	functions.read_time_stamp(filename_stamp)
			print (colored("\tA previous command generated results on: %s [%s -- Bacteriophage]" %(stamp, name), 'yellow'))
		else:		
			## Call phispy
			bacteriophage.ident_bacteriophage(gbk_file[0], name, outdir_phage, options.training_set, Debug, 
												window_size=options.window_size, number_phage_genes=options.phage_genes)
			
			## Parse results
			bacteriophage.results_PhiSpy(outdir_phage, name)
			
		## save results
		phage_Results[name] = outdir_phage
		
	####################
	## Genomic Island analysis
	####################
	if genomic_island_bool:

		#functions.print_sepLine("*",50, False)
		#print (' Genomic island identification analysis')
		#functions.print_sepLine("*",50, False)

		## for each sample
		outdir_GI = functions.create_subfolder("genomic_island", output_dir)

		## debug message
		if (Debug):
			print (colored("**DEBUG: Dir: "+ str(outdir_GI), 'yellow'))
		
		## Call phispy
		genomic_island.GI_module(gbk_file[0], name, outdir_GI, Debug, options.cutoff_dinuc_bias, options.min_length)
		
		## save results
		GI_Results[name] = outdir_GI
	
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

