#usr/bin/en python3
'''
This code prepares samples for further analysis.
Jose F. Sanchez
Copyright (C) 2019 Lauro Sumoy Lab, IGTP, Spain
'''
## useful imports
import time
import io
import os
import re
import sys

## import my modules
from BacterialTyper import functions
from BacterialTyper import config

###############
def select_samples (list_samples, samples_prefix, pair=True, exclude=False):
    
    #Get all files in the folder "path_to_samples"    
	sample_list = []
	for names in samples_prefix:
		for path_fastq in list_samples:	
			fastq = os.path.basename(path_fastq)
			samplename_search = re.search(r"(%s).*" % names, fastq)
			enter = ""
			if samplename_search:
				if (exclude): ## exclude==True
					enter = False
				else: ## exclude==True
					enter = True
			else:
				if (exclude): ## exclude==True
					enter = True
				else: ## exclude==True
					enter = False
					
			if enter:
				if fastq.endswith('.gz'):
					sample_list.append(path_fastq)
				elif fastq.endswith('fastq'):
					sample_list.append(path_fastq)
				else:
					print ("** ERROR: ", path_fastq, 'is a file that is neither in fastq.gz or .fastq format, so it is not included')

	## discard duplicates if any
	non_duplicate_samples = list(set(sample_list))	
	discard_samples = []
	
	if (pair):
		for path_files in non_duplicate_samples:
			files = os.path.basename(path_files)
			name_search = re.search(r"(.*)(\_1|_2)(\.f.*q)(\..*)", files)
			if name_search:
				file_name = name_search.group(1)
				read_pair = name_search.group(2)
				ext = name_search.group(3)
				gz = name_search.group(4)
				dirN = os.path.dirname(path_files)

				## get second pair
				paired = ""
				if (read_pair == '_1'):
					paired = dirN + '/' + file_name + '_2' + ext
				else:
					paired = dirN + '/' + file_name + '_1' + ext
	
				if gz:
					unzip = dirN + '/' + file_name + read_pair + ext
					if os.path.isfile(unzip):
						discard_samples.append(path_files)
						if os.path.isfile(paired):
							## both files are unzipped and available
							discard_samples.append(paired)
							print (file_name + '\t' + unzip + '\t&\t' + paired)							
							## save both files
					else:
						## gunzipped
						if os.path.isfile(paired + gz):
							## save both files
							print (file_name + '\t' + path_files + '\t&\t' + paired + gz)
				else:
					## not gunzipped files					
					if os.path.isfile(paired):
						## save both files
						print (file_name + '\t' + path_files + '\t&\t' + paired)
			
	exit()	
	for samples in discard_samples:
		non_duplicate_samples.remove(samples)
						
	non_duplicate_samples2 = list(set(non_duplicate_samples))						
	number_samples = len(non_duplicate_samples2)	

	print ("\t\t- ", number_samples," samples selected")
	return sorted(non_duplicate_samples)
	
###############

###############    
def one_file_per_sample(final_sample_list, path_to_samples, directory, read, output_file, prefix_list, num_threads):
	## merge sequencing files for sample, no matter of sector or lane generated.
	grouped_subsamples = []
	bigfile_list = []
	commands2sent = []
	
	output_file = open(output_file, 'a')
	output_file.write("\nMerge samples:\n")
	
	for samplex in final_sample_list:
		if samplex not in grouped_subsamples:
			samplename = []
			subsamples = []
			for prefix in prefix_list:
				samplename_search = re.search(r"(%s)\_(.*)" % prefix, samplex)
				if samplename_search:			
					original_name = samplename_search.group(1)
					commonname = original_name + "_" + read + ".fastq"
					bigfilepath = directory + "/" + commonname
					bigfile_list.append(commonname)	
					
					for sampley in final_sample_list:
						if original_name in sampley:
							subsamples.append(path_to_samples + "/" + sampley)
							grouped_subsamples.append(sampley)
					if not os.path.isfile(bigfilepath) or os.stat(bigfilepath).st_size == 0:
						partsofsample = ' '.join(sorted(subsamples))
						cmd = 'cat %s >> %s' %(partsofsample, bigfilepath)						
						## DUMP in file					
						output_file.write(cmd)   
						output_file.write('\n')						
						## get command				
						commands2sent.append(cmd)
					else:
						print ('\t + Sample %s is already merged' % commonname)
	
	## close file
	output_file.close()		
	#sent commands on threads
	functions.sender(commands2sent, num_threads)	
	print ('There are' , len(bigfile_list) , 'samples after merging for read' , read, '\n')
	return bigfile_list
###############

