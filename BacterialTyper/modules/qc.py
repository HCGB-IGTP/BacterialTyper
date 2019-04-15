#!/usr/bin/env python3
'''
This code calls 
Jose F. Sanchez
Copyright (C) 2019 Lauro Sumoy Lab, IGTP, Spain
'''

## import my modules
from BacterialTyper import trimmomatic_call
from BacterialTyper import multiQC_report
from BacterialTyper import spades_assembler

####
def fastqc(options):
	print ("Hi there! lets call FASTQC\n")	

def assembly_check(options):
	print ("Hi there! lets call BUSCO\n")

	# spades_assembler.stats(fasta_file)
