#!/usr/bin/env python3
##########################################################
## Jose F. Sanchez										##
## Copyright (C) 2019 Lauro Sumoy Lab, IGTP, Spain		##
##########################################################
"""
Generates a search among NCBI taxonomy and retrieves available genomes of interest
"""

## Original code extracted from: ncbi-genome-download contribution script
## https://raw.githubusercontent.com/kblin/ncbi-genome-download/master/contrib/gimme_taxa.py

## useful imports
import time
import io
import os
import re
import sys
import pandas as pd

from sys import argv
from io import open
from termcolor import colored
from ete3 import NCBITaxa

## import my modules
from BacterialTyper.scripts import functions
from BacterialTyper.config import set_config
from BacterialTyper.scripts import edirect_caller

##########################################################################################
def NCBItaxa_db(dbfile_path):
	
	# Generate the ete3 NCBI taxa object
	ncbi = NCBITaxa(dbfile_path)

	## check days passed
	#ncbi.update_taxonomy_database()
	
	
##########################################################################################
def get_descendant_taxid(taxid_name, dbfile_path):
	## use NCBI taxa class from ete3

	# Get NCBI taxa object
	ncbi = NCBITaxa_db(dbfile_path)
	
	## Look for taxid name provided
	taxid = ncbi.get_name_translator([taxid_name])[taxid_name][0]

	# Get all taxa from a given group
	descendent_taxa = ncbi.get_descendant_taxa(taxid)
	descendent_taxa_names = ncbi.translate_to_names(descendent_taxa)


