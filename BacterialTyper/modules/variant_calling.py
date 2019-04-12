 #!/usr/bin/env python3
'''
This code calls KMerGenie for k-mer optimization and Mccortex for the variant calling
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

## import my modules
from BacterialTyper import functions
from BacterialTyper import config

## https://samtools.github.io/hts-specs/VCFv4.2.pdf
## bcftools 

## /software/debian-8/bio/bcftools/bin/plot-vcfstats
## bcftools stats -s WTCHG_370809_205154 bubbles.joint.links.k51.geno.vcf > test
## plot-vcfstats test

## https://pyvcf.readthedocs.io/en/latest/

######
def	help_options():
	print ("\nUSAGE: python %s data reference name mmcortex_bin kmergenie_bin threads path\n"  %os.path.realpath(__file__))

######
def main():

  	## control if options provided or help
	if len(sys.argv) > 1:
		print ("")
	else:
		help_options()
		exit()    	

	data = os.path.abspath(argv[1])
	reference = os.path.abspath(argv[2])
	sample = argv[3]
	mccortex_bin = argv[4]
	kmergenie_bin = argv[5]
	threads = int(argv[6])
	path = 	argv[7]
	
	
	
######

'''******************************************'''
if __name__== "__main__":
	main()



