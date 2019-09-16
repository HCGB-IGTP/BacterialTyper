#usr/bin/env python
'''
This code calls mash software to... 
Jose F. Sanchez
Copyright (C) 2019 Lauro Sumoy Lab, IGTP, Spain
'''
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
import shutil

## import my modules
from BacterialTyper import functions
from BacterialTyper import config
from BacterialTyper import database_generator

##################################################
def helpMash():
	print (colored("\n\n***** TODO: Generate this help message *****\n\n", 'red'))

##################################################
def create_mash_database(option, name, database_folder, Debug):

	##############################################################
	## update Mash database with user_data information retrieved
	##############################################################
	mash_bin = config.get_exe('mash')
	
	print ('\n\n+ Update Mash database with information retrieved')
	print ('+ Get information...')

	option_string = ""
	if option == 'NCBI':
		option_string = 'genbank'
		
	## read database 
	db_frame = database_generator.getdbs(option, database_folder, option_string, Debug)
	data_db = database_generator.get_database(db_frame, Debug)
	
	## debug message
	if Debug:
		print (colored("**DEBUG: Retrieve database (%s)" %option + " **", 'yellow'))
		print (db_frame)
		
		print (colored("**DEBUG: Get database (%s)" %option + " **", 'yellow'))
		print (data_db)

	
	list_of_files = data_db['genome'].tolist()

	mash_folder = functions.create_subfolder('Mash_db', database_folder)	
	mash_db = functions.create_subfolder(name, mash_folder)	
	
	#####
	print ('+ Prepare database: ', name)
	## 
	lineList = []
	toIndexList = []
	indexedList = []		

	###
	## read db in fold_name and get index files
	info = mash_db + '/' + name + '.db'
	if os.path.exists(info):
		lineList = functions.readList_fromFile(info)
		option = 'add'
		
	for f in list_of_files:
		file_name = f
		baseName = os.path.basename(file_name)			
		## check if already index
		if baseName in lineList:
			print (colored('+ File %s is already available in database %s' %(baseName, name), 'green'))
			indexedList.append(file_name)
		else:
			toIndexList.append(file_name)		
	
	#########################################
	if toIndexList:

		## generate batch and call
		info2 = mash_db + '/.batch_entries.txt'
		functions.printList2file(info2, toIndexList)
		
		### call min hash
		## debug message
		if (Debug):
			print (colored("**DEBUG: Database (%s) is indexed" %name + " **", 'yellow'))
			print (mash_db)
			
		status = sketch_database(info2, mash_bin, name)
		
		## finish
		final_list = set(lineList + toIndexList + indexedList)
		final_list_name = [os.path.basename(f) for f in final_list]
		
		functions.printList2file(info, final_list_name)
		count_files = len(toIndexList)
		print ('+ %s samples have been added to the database' %count_files)

	else:
		print ('\n+ No new sequences were added to the database.')
		return (mash_db + '/' + name)			

	## check if previously indexed
	if (status): #true
		return (mash_db + '/' + name)			
	else: #false
		## debug message
		if (Debug):
			print (colored("**DEBUG: Database (%s) is not indexed" %name + " **", 'yellow'))
			return (False)
	

##################################################		
def sketch_database(list_files, mash_bin, out_file, name, folder):
	
	### -p threads
	mash_cmd = ""
	if len(list_files) == 1:
		print ('\t+ Skecthing sample: ', name)
		mash_cmd = '%s sketch -o %s %s' %(mash_bin, out_file, list_files[0])

	else:
		print ('\t+ Skecthing list of samples provided: ', name)
		## batch
		out_batch = folder + "/.batch_entries.txt"
		functions.printList2file(out_batch, list_files)
		mash_cmd = '%s sketch -l -o %s %s' %(mash_bin, out_file, out_batch)
	
	return(functions.system_call(mash_cmd))

##################################################		
def triangle_matrix_generation():
	print ()
	

##################################################		
def tree_construction():
	print ()


##################################################
def	help_options():
	print ("\nUSAGE: python %s\n"  %os.path.realpath(__file__))

##################################################
def main():
	## this code runs when call as a single script

  	## control if options provided or help
	if len(sys.argv) > 1:
		print ("")
	else:
		help_options()
		exit()
	
	## to do: implement main function
	mash_bin = config.get_exe('mash')
	
	print (mash_bin)
	

##################################################
if __name__== "__main__":
	main()
