#usr/bin/env python
'''
This code downloads information for assembly IDs provided
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
import pandas as pd
import ncbi_genome_download as ngd

######
def NCBIdownload(data2donwload):
	print (data2donwload[0]['NCBI_assembly_ID'])
	

######
def get_data(ID_file):	
	print ("+ Obtaining information from file:")
	print (ID_file)
	data = pd.read_csv(ID_file, header=0)

	print ("\n+ Data:")
	print (data)
	
	return(data)

######
def main():
	## this code runs when call as a single script

  	## control if options provided or help
	if len(sys.argv) > 1:
		print ("")
	else:
		help_options()
		exit()
	
	## get arguments provided
	ID_file = os.path.abspath(argv[1])
	
	## get file information
	strains2get = get_data(ID_file)
	
	## download
	NCBIdownload(strains2get)
	

######
def help_options():
	print ("\nUSAGE: python %s ID_file\n"  %os.path.realpath(__file__))
	print ("---------------")
	print ("ID_file example")
	print ("CSV file containing: genus,species name,strain name,Genbank or Refseq ID")	
	print ("---------------")
	print ("Genus,species,name,NCBI_assembly_ID")
	print ("genus,species,strain1,GCx_0000001")
	print ("genus,species,strain2,GCx_0000002")
	print ("genus,species,strain3,GCx_0000003")
	print ("genus,species,strain4,GCx_0000004")
	print ("genus,species,strain5,GCx_0000005")
	print ("*******")
	print ("______________")
		
######
'''******************************************'''
if __name__== "__main__":
	main()

