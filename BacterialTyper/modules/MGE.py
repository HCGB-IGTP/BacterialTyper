#!/usr/bin/env python3
'''
This code calls several programs to generate a Mobile Genetic Element analysis
Jose F. Sanchez
Copyright (C) 2019 Lauro Sumoy Lab, IGTP, Spain
'''
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
from BacterialTyper import functions
from BacterialTyper import config
from BacterialTyper import annotation
from BacterialTyper import bacteriophage
from BacterialTyper.modules import sample_prepare
import PhiSpy

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
	functions.boxymcboxface("Mobile Genetic Elements (MGE) identification")

	print ("--------- Starting Process ---------")
	functions.print_time()
	
	## print information and exit
	if options.help_PhiSpy:
		PhiSpy.print_list()
		exit()
		
	## absolute path out
	outdir = os.path.abspath(options.output_folder)

	## generate output folder
	functions.create_folder(outdir)

	### time stamp
	start_time_partial = start_time_total
	start_time_partial_assembly = start_time_partial
	
	###
	pd_samples_retrieved = pd.DataFrame()
	
	### batch provided
	if options.batch:
		## csv file containing sample name and file path
		pd_samples_retrieved = pd.read_csv(options.batch, sep=',',header=None)
		pd_samples_retrieved.columns = ["samples", "tag", "file"]
	
	else:
		fasta_ext = ('.fna', '.fasta')
		## assemblies
		assemblies_pd_samples_retrieved = pd.DataFrame()
		if options.input_assemblies:
			input_dir_assemblies = os.path.abspath(options.input_assemblies)

			print ("+ Retrieve all genomes assembled...")
			assemblies_pd_samples_retrieved = sample_prepare.get_files(options, input_dir_assemblies, "chromosome", fasta_ext)

		## plasmid
		plasmid_pd_samples_retrieved = pd.DataFrame()
		if options.input_plasmids:
			input_dir_plasmids = os.path.abspath(options.input_plasmids)

			print ("+ Retrieve all plasmids assembled...")
			plasmid_pd_samples_retrieved = sample_prepare.get_files(options, input_dir_plasmids, "plasmid", fasta_ext)
	
		## protein
		annotations_pd_samples_retrieved = pd.DataFrame()
		if options.input_annotations:
			input_dir_annotations = os.path.abspath(options.input_annotations)

			print ("+ Retrieve all annotations generated...")
			annotations_pd_samples_retrieved = sample_prepare.get_files(options, input_dir_annotations, "annotation", ('.gbf', '.faa', '.gff'))

		#print (assemblies_pd_samples_retrieved)
		#print (plasmid_pd_samples_retrieved)
		#print (annotations_pd_samples_retrieved)
		
		frames = [assemblies_pd_samples_retrieved, plasmid_pd_samples_retrieved, annotations_pd_samples_retrieved]
		pd_samples_retrieved = pd.concat(frames, ignore_index=True)

	print (pd_samples_retrieved)

	exit()

	

