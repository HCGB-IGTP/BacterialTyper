#!/usr/bin/env python3
##########################################################
## Jose F. Sanchez										##
## Copyright (C) 2019 Lauro Sumoy Lab, IGTP, Spain		##
##########################################################
"""
Prints an index of citation for the different packages and other softwares employed here.
"""

## import my modules
from BacterialTyper.config import set_config
import HCGB.functions.aesthetics_functions as HCGB_aes
from BacterialTyper import __version__ as pipeline_version

from termcolor import colored
import os, sys
import pandas as pd

#########
def software_citation():
	soft_citation ={
		#'CARD':['The Comprehensive Antibiotic Resistance Database', 'McArthur et al 2013', 'PMID: 23650175', 'https://card.mcmaster.ca/'],
		'ARIBA':['-', '-', '-', '-'],
		'AUGUSTUS':['-', '-', '-', '-'],
		'Bowtie2':['-', '-', '-', '-'],
		'BLAST':['-', '-', '-', '-'],
		'cdhit':['-', '-', '-', '-'],
		'FASTQC':['-', '-', '-', '-'],
		'KMA':['-', '-', '-', '-'],
		'PROKKA':['-', '-', '-', '-'],
		'SPADES':['-', '-', '-', '-'],
		'TRIMMOMATIC':['-', '-', '-', '-']		
	}
	return(soft_citation)

#########
def ariba_citation():
	ARIBA_DB_citation ={
		'CARD':['The Comprehensive Antibiotic Resistance Database', 'McArthur et al 2013', 'PMID: 23650175', 'https://card.mcmaster.ca/'],
		'MEGARes':['MEGARes: an antimicrobial database for high throughput sequencing', 'Lakin et al 2016', 'PMID: 27899569', 'http://megares.meglab.org/'],
		'PlasmidFinder':['PlasmidFinder and pMLST: in silico detection and typing of plasmids', 'Carattoli et al 2014', 'PMID: 24777092', '-'],
		'ResFinder':['Identification of acquired antimicrobial resistance genes', 'Zankari et al 2012', 'PMID: 22782487', '-'],
		'ARG-ANNOT':['ARG-ANNOT: a new bioinformatic tool to discover antibiotic resistance genes in bacterial genomes', 'Gupta et al 2014', 'PMID: 24145532', '-'],
		'srst2':['SRST2: Rapid genomic surveillance for public health and hospital microbiology labs', 'Inouye et al 2014', 'PMID: 25422674', 'https://github.com/katholt/srst2'],
		'VFDB':['VFDB 2019: a comparative pathogenomic platform with an interactive web interface.', 'Liu et al 2019', 'PMID: 30395255', 'http://www.mgc.ac.cn/VFs/main.htm'],
		'VirulenceFinder':['Real-time whole-genome sequencing for routine typing surveillance and outbreak detection of verotoxigenic Escherichia coli', 'Joensen al 2014', 'PMID: 24574290', '-']
	}
	return(ARIBA_DB_citation)

#########
def print_all():
	print ("")
	HCGB_aes.print_sepLine("+", 50, 'yellow')
	print ("\tSOFTWARE")
	HCGB_aes.print_sepLine("+", 50, 'yellow')
	print ("Third party softwares included or employed during the pipeline workflow.")
	print ("")
	df_software_citation = pd.DataFrame.from_dict(software_citation(), orient='index', columns=('Article Title', 'Authors', 'PUBMED ID', 'Website'))	
	df_software_citation.index.names = ['Software']
	pd.set_option('display.max_colwidth', None)
	pd.set_option('display.max_columns', None)
	print (df_software_citation)
	print ("")
	
	HCGB_aes.print_sepLine("+", 50, 'yellow')
	print ("\tDATABASES")
	HCGB_aes.print_sepLine("+", 50, 'yellow')
	print ("")
	print ("Please cite according to your selection.")
	print ("")
	
	HCGB_aes.print_sepLine("+", 50, False)
	print ("\tARIBA databases")
	HCGB_aes.print_sepLine("*", 50, False)	
	df_ARIBA_DB_citation = pd.DataFrame.from_dict(ariba_citation(), orient='index', columns=('Article Title', 'Authors', 'PUBMED ID', 'Website'))	
	df_ARIBA_DB_citation.index.names = ['Databases']
	print (df_ARIBA_DB_citation)
	print ("\n")
	
	HCGB_aes.print_sepLine("*", 50, False)
	print ("\tKMA software & databases")
	HCGB_aes.print_sepLine("*", 50, False)

	print ()	
	print ()	

	HCGB_aes.print_sepLine("*", 50, False)
	print ("\tBUSCO software & dataset")
	HCGB_aes.print_sepLine("*", 50, False)
	print ("BUSCO applications from quality assessments to gene prediction and phylogenomics.")
	print ("Robert M. Waterhouse, Mathieu Seppey, Felipe A. Simão, Mose Manni, Panagiotis ")
	print ("Ioannidis, Guennadi Klioutchnikov, Evgenia V. Kriventseva, and Evgeny M. Zdobnov")
	print ("Mol Biol Evol, published online Dec 6, 2017, doi: 10.1093/molbev/msx319")
	print ()	
	print ("BUSCO: assessing genome assembly and annotation completeness with single-copy orthologs.")
	print ("Felipe A. Simão, Robert M. Waterhouse, Panagiotis Ioannidis, Evgenia ")
	print ("V. Kriventseva, and Evgeny M. Zdobnov")
	print ("Bioinformatics, published online June 9, 2015, doi: 10.1093/bioinformatics/btv351")
	print ()	
	print ("For further details, please visit: https://busco.ezlab.org/ or https://www.orthodb.org/")

	print ()	
	print ()		
	
	
#########
def only_us():
	print ("Citation to be included when properly finished...\n\n")

#########
def run(options):
	HCGB_aes.pipeline_header("BacterialTyper", ver=pipeline_version)
	HCGB_aes.boxymcboxface("Citation")
	if (options.option == 'all'):
		print (colored("\n+ BacterialTyper citation:", 'blue'))
		only_us()

		print (colored("\n+ Other softwares employed in the pipeline", 'blue'))
		print_all()
		
	elif (options.option == 'only'):
		print (colored("\n+ BacterialTyper citation:", 'blue'))
		only_us()
	
	return()
	
