#!/usr/bin/env python3
'''
This code calls different submodules to assemble each sample using SPADES, check quality using BUSCO and annotate using PROKKA
Jose F. Sanchez
Copyright (C) 2019 Lauro Sumoy Lab, IGTP, Spain
'''

## import my modules
from BacterialTyper import trimmomatic_call
from BacterialTyper import spades_assembler
from BacterialTyper import config
from BacterialTyper import annotation
from BacterialTyper.modules import qc

####
def run(options):

	### 

	print ("Hi there! lets assemble data using SPADES\n")


	## qc.assembly_check()

####
##def assembly_check(options):
##	print ("Hi there! lets call BUSCO\n")
