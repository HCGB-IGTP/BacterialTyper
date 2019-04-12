#usr/bin/env python
'''
This code calls KMA software to find the best match (species identification) to the reads in one
or more fastq files or one fasta file in a (kmer) database produced using the KMA program (Philip T.L.C. Clausen, 
Frank M. Aarestrup & Ole Lund, "Rapid and precise alignment of raw reads against redundant databases with KMA", 
BMC Bioinformatics, 2018;19:307.) 
Jose F. Sanchez
Copyright (C) 2019 Lauro Sumoy Lab, IGTP, Spain
'''
## useful imports
import time
import io
import os
import re
import sys
from sys import argv
from io import open

## import modules
pythonDir = os.path.dirname(os.path.realpath(__file__)) + '/../tools/python'
sys.path.append(pythonDir)
import functions

## config
configDir = os.path.dirname(os.path.realpath(__file__)) + '/../config/'
sys.path.append(configDir)
import config

##########
def index_database(database_entries, kma_bin, index_name, option):
	
	########################################################################################
	## 								KMA_index-1.2.2							
	########################################################################################
	# kma_index creates the databases needed to run KMA, from a list of fasta files given.
	# Options are:		
	#				Desc:									Default:
	#
	#	-i			Input/query file name (STDIN: "--")		None
	#	-o			Output file								Input/template file
	#	-batch		Batch input file
	#	-deCon		File with contamination (STDIN: "--")	None/False
	#	-batchD		Batch decon file
	#	-t_db		Add to existing DB						None/False
	#	-k			Kmersize								16
	#	-k_t		Kmersize for template identification	16
	#	-k_i		Kmersize for indexing					16
	#	-ML			Minimum length of templates	kmersize 	(16)	
	#	-CS			Start Chain size						1 M
	#	-ME			Mega DB									False
	#	-NI			Do not dump *.index.b					False
	#	-Sparse		Make Sparse DB ('-' for no prefix)		None/False
	#	-ht			Homology template						1.0
	#	-hq			Homology query							1.0
	#	-and		Both homolgy thresholds
	#				has to be reached						or
	#	-v			Version
	#	-h			Shows this help message
	#######################################################################################
	
	if (option == "new"):
		cmd_kma_index = "%s -batch %s -o %s" %(kma_bin, database_entries, index_name)
	elif (option == "add"):
		cmd_kma_index = "%s -batch %s -o %s -t_db" %(kma_bin, database_entries, index_name)

##########	
def kmerfinder_call():
	print()
	
##########
def kmerfinder_module():
	print()

######
def	help_options():
	print ("\nUSAGE: python %s file1 file2 name SPADES_bin threads path\n"  %os.path.realpath(__file__))

######
def main():
	## this code runs when call as a single script

  	## control if options provided or help
	if len(sys.argv) > 1:
		print ("")
	else:
		help_options()
		exit()    	


######

'''******************************************'''
if __name__== "__main__":
	main()
	
	
################################################################################################
## 										KMA v1.2.2
################################################################################################	
# KMA-1.2.2 mapps raw reads to a template database.
# Options are:		
##				Desc:													Default:	Requirements:
#
#	-o			Output file												None		REQUIRED
#	-t_db		Template DB												None		REQUIRED
#	-i			Input file name(s)										STDIN
#	-ipe		Input paired end file name(s)
#	-int		Input interleaved file name(s)
#	-k			Kmersize												DB defined
#	-e			evalue													0.05
#	-ConClave	ConClave version										1
#	-mem_mode	Use kmers to choose best template, and save memory		False
#	-ex_mode	Searh kmers exhaustively								False
#	-ef			Print additional features								False
#	-vcf		Make vcf file, 2 to apply FT							False/0
#	-deCon		Remove contamination									False
#	-dense		Do not allow insertions in assembly						False
#	-ref_fsa	Consensus sequnce will have "n" instead of gaps			False
#	-matrix		Print assembly matrix									False
#	-a			Print all best mappings									False
#	-mp			Minimum phred score										20
#	-5p			Cut a constant number of nucleotides from the 5 prime.	0
#	-Sparse		Only count kmers										False
#	-Mt1		Map only to "num" template.								0 / False
#	-ID			Minimum ID												1.0%
#	-ss			Sparse sorting (q,c,d)									q
#	-pm			Pairing method (p,u,f)									u
#	-fpm		Fine Pairing method (p,u,f)								u
#	-apm		Sets both pm and fpm									u
#	-shm		Use shared DB made by kma_shm							0 (lvl)
#	-1t1		Skip HMM												False
#	-ck			Count kmers instead of pseudo alignment					False
#	-ca			Make circular alignments								False
#	-boot		Bootstrap sequence										False
#	-bc			Base calls should be significantly overrepresented.		[True]
#	-bc90		Base calls should be both  significantly 
#				overrepresented, and have 90% agreement.				False
#	-bcNano		Call bases at suspicious deletions, made for nanopore.	False
#	-bcd		Minimum depth at base									1
#	-bcg		Maintain insignificant gaps
#	-and		Both mrs and p_value thresholds has 
#				to reached to in order to report a template hit
#	-mq			Minimum mapping quality									0
#	-mrs		Minimum alignment score, normalized to 
#				alignment length										0.50
#	-reward		Score for match											1
#	-penalty	Penalty for mismatch									-2
#	-gapopen	Penalty for gap opening									-3
#	-gapextend	Penalty for gap extension								-1
#	-per		Reward for pairing reads								7
#	-cge		Set CGE penalties and rewards							False
#	-t			Number of threads										1
#	-v			Version
#	-h			Shows this help message
#################################################################################################

## kma index -batch ../database.txt -k 30 -k_t 30 -k_i 30 -o index_db_k30
## kma -ipe WTCHG_370809_205154/WTCHG_370809_205154_trim_R1.fastq WTCHG_370809_205154/WTCHG_370809_205154_trim_R2.fastq -o kma_search3 -t_db index_kma/index_db_k30 -t 3 -apm u -1t1 -k 30
	

