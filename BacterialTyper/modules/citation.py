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

#########
def print_all():

	functions.print_sepLine("+", 50)
	print ("ARIBA databases")
	functions.print_sepLine("+", 50)
	print ("card -> 'The Comprehensive Antibiotic Resistance Database', McArthur et al 2013, PMID: 23650175")
	print ("megares  ->  MEGARes: an antimicrobial database for high throughput sequencing', Lakin et al 2016, PMID: PMC5210519")
	print ("plasmidfinder  ->  'PlasmidFinder and pMLST: in silico detection and typing of plasmids', Carattoli et al 2014, PMID: 24777092")
	print ("resfinder  ->  'Identification of acquired antimicrobial resistance genes', Zankari et al 2012, PMID: 22782487")
	print ("srst2_argannot  ->  'ARG-ANNOT, a new bioinformatic tool to discover antibiotic resistance genes in bacterial genomes', Gupta et al 2014, PMID: 24145532")
	print ("vfdb_core & vfdb_full  ->  'VFDB 2016: hierarchical and refined dataset for big data analysis-10 years on',Chen LH et al 2016, Nucleic Acids Res. 44(Database issue):D694-D697. PMID: 26578559")
	print ("virulencefinder  ->  'Real-time whole-genome sequencing for routine typing, surveillance, and outbreak detection of verotoxigenic Escherichia coli', Joensen al 2014, PMID: 24574290")
	print ("\n")
	
#########
def only_us():
	print ("Citation to be included when properly finished...\n\n")

#########
def run(options):
	if (options.all):
		print (colored("+ Other softwares employed in the pipeline", 'blue'))
		print_all()
		
		print (colored("\n+ BacterialTyper citation:", 'blue'))
		only_us()
		
	else:
		only_us()
	
	
