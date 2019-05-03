#!/usr/bin/env python3
'''
This code generates prints an index of citation for the different packages and other softwares employed here.
Jose F. Sanchez
Copyright (C) 2019 Lauro Sumoy Lab, IGTP, Spain
'''

## import my modules
from BacterialTyper import config, functions
from termcolor import colored
import os, sys
import pandas as pd

#########
def ariba_citation():
	ARIBA_DB_citation ={
		'CARD':['The Comprehensive Antibiotic Resistance Database', 'McArthur et al 2013', 'PMID: 23650175'],
		'MEGARes':['MEGARes: an antimicrobial database for high throughput sequencing', 'Lakin et al 2016', 'PMID: 27899569'],
		'PlasmidFinder':['PlasmidFinder and pMLST: in silico detection and typing of plasmids', 'Carattoli et al 2014', 'PMID: 24777092'],
		'ResFinder':['Identification of acquired antimicrobial resistance genes', 'Zankari et al 2012', 'PMID: 22782487'],
		'ARG-ANNOT':['ARG-ANNOT: a new bioinformatic tool to discover antibiotic resistance genes in bacterial genomes', 'Gupta et al 2014', 'PMID: 24145532'],
		'srst2':['SRST2: Rapid genomic surveillance for public health and hospital microbiology labs', 'Inouye et al 2014', 'PMID: 25422674'],
		'VFDB':['VFDB 2019: a comparative pathogenomic platform with an interactive web interface.', 'Liu et al 2019', 'PMID: 30395255'],
		'VirulenceFinder':['Real-time whole-genome sequencing for routine typing surveillance and outbreak detection of verotoxigenic Escherichia coli', 'Joensen al 2014', 'PMID: 24574290']
	}
	return(ARIBA_DB_citation)

#########
def print_all():
	functions.pipeline_header()
	functions.print_sepLine("+", 50)
	print ("\tARIBA databases")
	functions.print_sepLine("+", 50)
	
	print ("- Please cite according to your selection.")
	
	df_ARIBA_DB_citation = pd.DataFrame.from_dict(ariba_citation(), orient='index')	
	df_ARIBA_DB_citation.index.names = ['Databases']
	
	pd.set_option('display.max_colwidth', -1)
	pd.set_option('display.max_columns', None)
	print (df_ARIBA_DB_citation)
	print ("\n")
	
#########
def only_us():
	print ("Citation to be included when properly finished...\n\n")

#########
def run(options):

	functions.boxymcboxface("Citation")
	if (options.option == 'all'):
		print (colored("\n+ BacterialTyper citation:", 'blue'))
		only_us()

		print (colored("\n+ Other softwares employed in the pipeline", 'blue'))
		print_all()
		
	elif (options.option == 'only'):
		print (colored("\n+ BacterialTyper citation:", 'blue'))
		only_us()
	
	
