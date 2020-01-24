#!/usr/bin/env python3
##########################################################
## Jose F. Sanchez										##
## Copyright (C) 2019 Lauro Sumoy Lab, IGTP, Spain		##
##########################################################
from _ast import If
from sphinx.search import tr
'''
Prepares samples for further analysis.
'''
## useful imports
import time
import io
import os
import re
import sys
import pandas as pd
from termcolor import colored
import shutil
import concurrent.futures

## import my modules
from BacterialTyper import functions
from BacterialTyper import config

###############
def help_format():

	functions.boxymcboxface("Name format for samples")

	print ("Format for fastq files can be:\n")
	functions.print_sepLine("*",20,"red")
	print ("[1] Single end files")
	functions.print_sepLine("*",20,"red")
	print ("name.fastq.gz")
	print ("name.fastq")
	print ("name.fq")
	print ("\n")

	functions.print_sepLine("*",20,"red")
	print ("[2] Paired-end files")
	functions.print_sepLine("*",20,"red")
	print ("Read1 => name_1.fastq.gz\tname_R1.fastq.gz")
	print ("Read2 => name_2.fastq.gz\tname_R2.fastq.gz")
	print ("\n")

	functions.print_sepLine("*",55,"red")
	print ("[3] Containing Lane information (*L00x* and/or *00x*):")
	functions.print_sepLine("*",55,"red")
	print ("name_L00x_R1.fastq.gz\tname_L00x_R2.fastq.gz")
	print ("name_L00x_1.fastq.gz\tname_L00x_2.fastq.gz")
	print ("name_L00x_R1_00x.fastq.gz\tname_L00x_R1_00x.fastq.gz")
	print (" ** Merge lane files using module option merge module prep **")
	print ("\n")

	functions.print_sepLine("*",15,"red")
	print ("[4] Extensions:")
	functions.print_sepLine("*",15,"red")
	print ("name_L00x_R2.fastq\tname_L00x_R2.fq\nname_L00x_R2.fastq.gz\tname_L00x_R2.fq.gz")
	print ("\n")

###############
def get_fields(file_name_list, pair=True, Debug=False):

	## init dataframe
	name_columns = ("sample", "dirname", "name", "lane", "read_pair","lane_file","ext","gz", "tag")
	name_frame = pd.DataFrame(columns=name_columns)
	
	## loop through list
	for path_files in file_name_list:
		
		## get file name
		file_name = os.path.basename(path_files)
		dirN = os.path.dirname(path_files)

		##	
		trim_search = re.search(r".*trim.*", file_name)
		lane_search = re.search(r".*\_L\d+\_.*", file_name)

		## get name
		if (pair):
		
			## pair could be: R1|R2 or 1|2
			## lane should contain L00x			
	
			if (trim_search):
				name_search = re.search(r"(.*)\_trim\_(R1|1|R2|2)\.(f.*q)(\..*){0,1}", file_name)
			else:
				## Lane files: need to merge by file_name: 33i_S5_L004_R1_001.fastq.gz
				if (lane_search):
					name_search = re.search(r"(.*)\_(L\d+)\_(R1|1|R2|2)(.*)\.(f.*q)(\..*){0,1}", file_name)
				else:
					name_search = re.search(r"(.*)\_(R1|1|R2|2)\.(f.*q)(\..*){0,1}", file_name)
		else:
			name_search = re.search(r"(.*)\.(f.*q)(\..*){0,1}", file_name)
	
		### declare
		name= ""
		lane_id= ""
		read_pair= ""
		lane_file= ""
		ext= ""
		gz= ""
	
		if name_search:
			name = name_search.group(1)
			if (pair):
				## Lane files: need to merge by file_name: 33i_S5_L004_R1_001.fastq.gz
				if (lane_search):
					lane_id = name_search.group(2)
					read_pair = name_search.group(3)
					lane_file = name_search.group(4)
					ext = name_search.group(5)
					gz = name_search.group(6)
				else:
					## could exist or not
					read_pair = name_search.group(2)
					ext = name_search.group(3)
					gz = name_search.group(4)
			else:
				ext = name_search.group(2)
				gz = name_search.group(3)
	
			name_frame.loc [len(name_frame)] = (path_files, dirN, name, lane_id, read_pair, lane_file, ext, gz, "reads")
	
		else:
			## debug message
			if (Debug):
				print (colored("**DEBUG: sampleParser.get_fields **", 'yellow'))
				print (colored("*** ATTENTION: Sample did not match the possible parsing options...", 'yellow'))
				print (file_name)

	return (name_frame)

###############
def select_samples (list_samples, samples_prefix, pair=True, exclude=False, Debug=False):
    
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
					## debug message
					if (Debug):
						print (colored("**DEBUG: sampleParser.select_samples **", 'yellow'))
						print (colored("** ERROR: %s is a file that is neither in fastq.gz or .fastq format, so it is not included" %path_fastq, 'yellow'))

	## discard duplicates if any
	non_duplicate_samples = list(set(sample_list))	
	
	## get fields
	name_frame_samples = get_fields(non_duplicate_samples, pair, Debug)	
	number_samples = name_frame_samples.index.size
	total_samples = set(name_frame_samples['name'].to_list())
	
	### get some stats
	if (number_samples == 0):
		print (colored("\n**ERROR: No samples were retrieved. Check the input provided\n",'red'))
		exit()
	print (colored("\t" + str(number_samples) + " files selected...", 'yellow'))
	print (colored("\t" + str(len(total_samples)) + " samples selected...", 'yellow'))
	if (pair):
		print (colored("\tPaired-end mode selected...", 'yellow'))
	else:
		print (colored("\tSingle end mode selected...", 'yellow'))
	
	## return info
	return (name_frame_samples)

###############
def select_other_samples (project, list_samples, samples_prefix, mode, extensions, exclude=False, Debug=False):

	## init dataframe
	name_columns = ("sample", "dirname", "name", "ext", "tag")

	## initiate dataframe
	df_samples = pd.DataFrame(columns=name_columns)

    #Get all files in the folder "path_to_samples"    
	sample_list = []
	for names in samples_prefix:
		for path_file in list_samples:	
			f = os.path.basename(path_file)
			dirN = os.path.dirname(path_file)
			#samplename_search = re.search(r"(%s).*" % names, f)
			samplename_search = re.search(r"(%s).*" % names, path_file)
			
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
				
				## project mode:
				if project:
					if mode == 'annot':
						#### /path/to/folder/annot/name.faa
						for ext in extensions:
							f_search = re.search(r".*\/%s\/(.*)\.%s$" %(mode, ext), path_file)
							if f_search:
								file_name = f_search.group(1) 
								df_samples.loc[len(df_samples)] = [path_file, dirN, file_name, ext, mode]	

					elif mode== 'assembly':
						#### name_assembly.faa
						f_search = re.search(r"(.*)\_%s\.%s$" %(mode, extensions), f)
						if f_search:
							file_name = f_search.group(1) 
							df_samples.loc[len(df_samples)] = [path_file, dirN, file_name, extensions, mode]	

					elif mode== 'mash':
						#### name.sig
						f_search = re.search(r".*\/%s\/(.*)\.%s$" %(mode, extensions[0]), path_file)
						if f_search:
							file_name = f_search.group(1) 
							df_samples.loc[len(df_samples)] = [path_file, dirN, file_name, extensions[0], mode]	

					else:
						f_search = re.search(r".*\/(.*)\/%s\/(.*)\_summary\.%s$" %(mode, extensions[0]), path_file)
						if f_search:
							### get information
							if mode == 'profile':
								name = f_search.group(1)
								db_name = f_search.group(2).split('_')[-1]
								if not name.startswith('report'):
									df_samples.loc[len(df_samples)] = [path_file, dirN, name, db_name, mode]	

							elif mode == 'ident':
								name = f_search.group(1)
								df_samples.loc[len(df_samples)] = [path_file, dirN, name, 'csv', mode]	

				## detached mode
				else:
					if f.endswith(extensions):
						file_name, ext = os.path.splitext(f)
						df_samples.loc[len(df_samples)] = [path_file, dirN, file_name, db_name, mode]	
					
	## debug message
	if (Debug):
		print (colored("**DEBUG: df_samples **", 'yellow'))
		print (df_samples)
	
	##
	number_samples = df_samples.index.size
	if (number_samples == 0):
		print (colored("\n**ERROR: No samples were retrieved for this option. Continue processing...\n",'red'))
		return (df_samples)
	print (colored("\t" + str(number_samples) + " samples selected from the input provided...", 'yellow'))

	return (df_samples)
	

###############    
def gunzip_merge(outfile, list_files):
	
	list_files = list(list_files)
	list_files.sort()
	print ("\tMerging files into: ", outfile)
	#print ("\tFiles: ", list_files)

	with open(outfile, 'wb') as wfp:
		for fn in list_files:
			with open(fn, 'rb') as rfp:
				shutil.copyfileobj(rfp, wfp)

	return()
	
###############    
def one_file_per_sample(dataFrame, outdir_dict, threads, outdir, Debug=False):
	## merge sequencing files for sample, no matter of sector or lane generated.
	
	list_samples = set(dataFrame['name'].tolist())
	print (colored("\t" + str(len(list_samples)) + " samples to be merged from the input provided...", 'yellow'))
	print ("+ Merging sequencing files for samples")

	##
	sample_frame = dataFrame.groupby(["name", "read_pair"])
	
	### get extension for files
	ext_list = dataFrame.ext.unique()
	gz_list = dataFrame.gz.unique()
	ext = ext_list[0] + gz_list[0] ## might generate a bug if several extension or some zip/unzip files provided

	# We can use a with statement to ensure threads are cleaned up promptly
	with concurrent.futures.ThreadPoolExecutor(max_workers=threads) as executor: ## need to do 1 by one as there is a problem with the working directory
		commandsSent = { executor.submit(gunzip_merge, outdir_dict[name[0]] + '/' + name[0] + '_' + name[1] + ext, sorted(set(cluster["sample"].tolist()))): name for name, cluster in sample_frame }
		for cmd2 in concurrent.futures.as_completed(commandsSent):
			details = commandsSent[cmd2]
			try:
				data = cmd2.result()
			except Exception as exc:
				print ('***ERROR:')
				print (cmd2)
				print('%r generated an exception: %s' % (details, exc))
				
							
	## return output name merged generated in dataframe
	name_columns = ("name", "dirname", "read_pair", "sample", "ext", "gz")
	name_frame = pd.DataFrame(columns=name_columns)

	## print to a file
	timestamp = functions.create_human_timestamp()
	merge_details = outdir + '/' + timestamp + '_prep_mergeDetails.txt'
	merge_details_hd = open(merge_details, 'w')

	for name, cluster in sample_frame: ## loop over samples
		outfile = outdir_dict[name[0]] + '/' + name[0] + '_' + name[1] + ext
		
		merge_details_hd.write("####################\n")		
		merge_details_hd.write("Sample: " + name[0] + '\n')
		merge_details_hd.write("Read: " + name[1] + '\n')
		merge_details_hd.write("Files:\n")
		merge_details_hd.write(",".join(cluster["sample"].tolist()))
		merge_details_hd.write('\n')
		merge_details_hd.write("####################\n")		
		
		name_frame.loc[len(name_frame)] = (name[0], outdir_dict[name[0]], name[1], outfile, ext_list[0], gz_list[0])

	merge_details_hd.close()
	return(name_frame)
	
###############
def help_options():
	print ("\nUSAGE:\npython %s in_folder out_folder threads\n"  %os.path.abspath(sys.argv[0]))

################################
def get_files(options, input_dir, mode, extension):
	
	## get list of input files
	files = []
	
	print ()
	functions.print_sepLine("-",50, False)
	print ('+ Getting files from input folder... ')
	print ('+ Mode: ', mode,'. Extension:', extension)
	if (options.project):
		### a folder containing a project is provided
		if os.path.exists(input_dir):
			#print ('+ Input folder exists')
			## get files in folder
			files = []
			for ext in extension:
				if mode == 'trim':
					files_tmp = functions.get_fullpath_list(input_dir)
					files = [s for s in files_tmp if ext in s]
				else:
					files_tmp = functions.retrieve_matching_files(input_dir, ext)				
					files.extend(files_tmp)
			
			files = set(files)
			
		else:
			## input folder does not exist...
			if (options.debug):
				print (colored("\n**DEBUG: sample_prepare.get_files input folder does not exists **", 'yellow'))
				print (input_dir)
				print ("\n")
			
			print (colored('***ERROR: input folder does not exist or it is not readable', 'red'))
			exit()
						
	else:
		### provide a single folder or a file with multiple paths (option batch)
		if (options.batch):
			if os.path.isfile(input_dir):
				dir_list = [line.rstrip('\n') for line in open(input_dir)]
				for d in dir_list:
					if os.path.exists(d):
						print ('+ Folder (%s) exists' %d)
						files = files + functions.get_fullpath_list(d)
					else:
						## input folder does not exist...
						if (options.debug):
							print (colored("\n**DEBUG: sample_prepare.get_files batch option; input folder does not exists **", 'yellow'))
							print (d)
							print ("\n")
						
		else:
			if os.path.exists(input_dir):
				print ('+ Input folder exists')
				## get files in folder
				files = functions.get_fullpath_list(input_dir)
			else:
				## input folder does not exist...
				if (options.debug):
					print (colored("\n**DEBUG: sample_prepare.get_files input folder does not exists **", 'yellow'))
					print (input_dir)
					print ("\n")
	
				print (colored('***ERROR: input folder does not exist or it is not readable', 'red'))
				exit()

	## get list of samples
	samples_names = []
	exclude=False

	if (options.in_sample):
		in_file = os.path.abspath(options.in_sample)
		samples_names = [line.rstrip('\n') for line in open(in_file)]
		print ('+ Retrieve selected samples to obtain from the list files available.')		
		exclude=False

		## in sample list...
		if (options.debug):
			print (colored("\n**DEBUG: sample_prepare.get_files include sample list **", 'yellow'))
			print (samples_names, '\n')


	elif (options.ex_sample):
		ex_file = os.path.abspath(options.ex_sample)
		samples_names = [line.rstrip('\n') for line in open(ex_file)]
		print ('+ Retrieve selected samples to exclude from the list files available.')		
		exclude=True

		## in sample list...
		if (options.debug):
			print (colored("\n**DEBUG: sample_prepare.get_files exclude sample list **", 'yellow'))
			print (samples_names, '\n')

	else:
		samples_names = ['.*']

	## discard some files obtain
	files = [s for s in files if 'single_copy_busco_sequences' not in s]
	files = [s for s in files if 'orphan' not in s]
	files = [s for s in files if 'augustus_output' not in s]
	files = [s for s in files if 'hmmer_output' not in s]
	files = [s for s in files if 'configs' not in s]
	files = [s for s in files if '00.0_0.cor.fastq.gz' not in s]
	files = [s for s in files if 'report_summary' not in s]
		
	## files list...
	if (options.debug):
		print (colored("\n**DEBUG: sample_prepare.get_files files list to check **", 'yellow'))
		##print ('DO NOT PRINT THIS LIST: It could be very large...')
		print (files, '\n')

	## get information
	if mode in ['fastq', 'trim']:
		pd_samples_retrieved = select_samples(files, samples_names, options.pair, exclude, options.debug)
	else:
		pd_samples_retrieved = select_other_samples(options.project, files, samples_names, mode, extension, exclude, options.debug)		
		
	return(pd_samples_retrieved)

######
def main():
	## this code runs when call as a single script

  	## control if options provided or help
	if len(sys.argv) > 1:
		print ("")
	else:
		help_options()
		exit()
	
	in_folder = os.path.abspath(sys.argv[1])
	out_folder = os.path.abspath(sys.argv[2])
	threads = int(sys.argv[3])
	
	## get files path
	list_files = []
	
	for root, dirs, files in os.walk(in_folder):
		for f in files:
			list_files.append(os.path.join(root,f))
	
	names = ['.*']
	datareturn = select_samples(list_files, names, True, False, True, True)
	one_file_per_sample(datareturn, out_folder, threads)
	
######
if __name__== "__main__":
	main()

