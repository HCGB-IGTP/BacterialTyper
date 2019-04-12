#usr/bin/env python
'''
This code calls fastqc analysis and parses results generated.
Jose F. Sanchez
Copyright (C) 2019 Lauro Sumoy Lab, IGTP, Spain
'''
## useful imports
import time
import io
import os
import re
import sys
import fastqcparser
from sys import argv
from io import open
import pandas as pd

import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
from pandas.plotting import table
from matplotlib.backends.backend_pdf import PdfPages

import numpy as np

## import my modules
from BacterialTyper import functions
from BacterialTyper import config

##################
class fastqcObject: 

    def __init__(self, sample=None, path=None, StatisticsFrame=None, StatusFrame=None):
        self.__sample=sample
        self.__path = path
        self.__StatusFrame = StatusFrame
        self.__StatisticsFrame = StatisticsFrame
        
    def print_header(self):
        if self.__sample:
            print("Working on sample: %s\nPath: %s " %(self.__sample, self.__path))
            
    def set_sample(self, sample):
        self.__sample = sample
        
    def get_sample(self):
        return self.__sample    

    def set_StatusFrame(self, by):
        self.__StatusFrame = by
        
    def get_StatusFrame(self):
        return self.__StatusFrame
    
    def set_StatisticsFrame(self, by):
        self.__StatisticsFrame = by
        
    def get_StatisticsFrame(self):
        return self.__StatisticsFrame
    
    def __repr__(self):
        return "fastqcObject('" + self.__path + "', Name: " +  str(self.__sample) + ")"

    def __str__(self):
        return "fastqcObject: " + self.__path + ", Name: " +  str(self.__sample)

############
def	help_options():
	print ("\nUSAGE:\npython %s folder file1 file2 name fastqc\n"  %os.path.realpath(__file__))

############
def call_fastqc(path, file1, file2, sample, fastqc_bin):	
	## call system for fastqc sample given
	name = functions.create_subfolder(sample, path)
	logFile = path + '/' + sample + '.log'
	
	if os.path.isfile(logFile):
		return ('OK')
	
	print ("+ Calling fastqc for samples...")	
	cmd_fastqc = '%s --extract -o %s %s %s > %s 2> %s' %(fastqc_bin, name, file1, file2, logFile, logFile)
	## send command	
	return (functions.system_call( cmd_fastqc ))
		
############
def parse_fastqcFile(resultsfile, name):
	## parse fastqc_data.txt file 
	print ("+ Parsing results from fastqc analysis")
	# load file
	f = fastqcparser.FastQCParser(resultsfile)
	
	# statistics
	col_names_statistics = ['filename', 'Sequences', 'Filtered', '%GC']
	statistics_df = pd.DataFrame(columns = col_names_statistics)
	
	sampleRead_search = re.search(r"(.*)_(R\d{1})(.*)", f.filename)
	read_pair = sampleRead_search.group(2)
	
	statistics_df.loc[0] = [f.filename, f.total_sequences,f.filtered_sequences,f.modules['Basic Statistics']['data'][-1][1]]
	#statistics_df.append([name, f.filename, f.encoding,f.total_sequences,f.filtered_sequences,f.modules['Basic Statistics']['data'][-1][1]], ignore_index=True)
	
	# status
	col_names_status = ['name', 'filename', 'readpair']
	modules_list = list(f.modules.keys())
	col_names_status.extend(modules_list)
	status_df = pd.DataFrame(columns = col_names_status)
	status = []
	
	for values in f.modules.keys():
		status.append(f.modules[values]['status'])
	
	status.insert(0, read_pair)
	status.insert(0, f.filename)
	status.insert(0, name)
	status_df.loc[0] = pd.Series(status).values

	return (statistics_df, status_df)
	
############
def run_module_fastqc(path, file1, file2, sample):	
	## Arguments provided via ARGVs
	fastqc_bin = config.CONFIGURATION['fastqc']
	codeReturn = call_fastqc(path, file1, file2, sample, fastqc_bin)
	if codeReturn == 'FAIL':
		exit()
	path_to_sample = path + '/' + sample
	return path_to_sample

############
def generateTable(dataFrame, fileName):
	## set index
	index_dataFrame = dataFrame.set_index(['name', 'readpair'])
	index_dataFrame.to_csv(fileName + '.txt', sep='\t')

	## plot table
	pdf_name = fileName + '.pdf'
	pp = PdfPages(pdf_name)

	fig,ax = plt.subplots()
	fig.patch.set_visible(False)
	
	ax.axis('off')
	ax.axis('tight')
	
	## color according to pass|warn|fail
	colors = index_dataFrame.applymap(lambda x: 'green' if x== 'pass' else ('yellow' if x== 'warn' else ('red' if x=='fail' else 'white' )))
	## https://www.rapidtables.com/web/color/html-color-codes.html

	tab = ax.table(
		cellText=index_dataFrame.values,
		#rowLabels=index_dataFrame.index.get_level_values(1), 
		rowLabels=index_dataFrame.index.values, 
		colLabels=index_dataFrame.columns,
		cellColours=colors.values,
		colWidths=[0.06 for x in index_dataFrame.columns],
		loc='center', colLoc = 'center', rowLoc='left', cellLoc='center')
	tab.auto_set_font_size(False)
	tab.set_fontsize(5)

	for cell in tab._cells:
		if cell[0]==0:
			tab._cells[cell].get_text().set_rotation(80)
			tab._cells[cell].set_height(0.25)
			tab._cells[cell].set_linewidth(0)
			tab._cells[cell]._loc = 'center'
			#tab._cells[cell].set_text_props(fontproperties=FontProperties(weight='bold'))
		
	fig.tight_layout()	
	pp.savefig()
	plt.close()
	pp.close()

def generateStatisticsTable(dataFrame, fileName):
	## set index
	index_dataFrame = dataFrame.set_index(['name', 'readpair'])
	index_dataFrame.to_csv(fileName + '.txt', sep='\t')

	## plot table
	pdf_name = fileName + '.pdf'
	pp = PdfPages(pdf_name)

	fig,ax = plt.subplots()
	fig.patch.set_visible(False)
	ax.axis('off')
	ax.axis('tight')
	
	tab = ax.table(
		cellText=index_dataFrame.values,
		#rowLabels=index_dataFrame.index.get_level_values(1), 
		rowLabels=index_dataFrame.index.values, 
		colLabels=index_dataFrame.columns,
		colWidths=[0.06 for x in index_dataFrame.columns],
		loc='center', colLoc = 'center', rowLoc='left', cellLoc='center')
	tab.auto_set_font_size(False)
	tab.set_fontsize(5)

	for cell in tab._cells:
		if cell[0]==0:
			tab._cells[cell].get_text().set_rotation(80)
			tab._cells[cell].set_height(0.25)
			tab._cells[cell].set_linewidth(0)
			tab._cells[cell]._loc = 'center'
			#tab._cells[cell].set_text_props(fontproperties=FontProperties(weight='bold'))
		
	fig.tight_layout()	
	pp.savefig()
	plt.close()
	pp.close()

############
def get_files(path):
	files = os.listdir(path)
	fastqc_files = []
	## get files generated with summary information
	for f in files:
		if f.endswith('_fastqc'):
			tmp = path + '/' + f + '/fastqc_data.txt'
			if os.path.isfile(tmp):
				fastqc_files.append(tmp)
	return (fastqc_files)

############
def main():
	## this code runs when call as a single script

  	## control if options provided or help
	if len(sys.argv) > 1:
		print ("")
	else:
		help_options()
		exit()    	

	path = os.path.abspath(argv[1])
	file1 = os.path.abspath(argv[2])
	file2 = os.path.abspath(argv[3])
	sample = argv[4]
	fastqc_bin = argv[5]

	##
	path_to_sample = run_module_fastqc(path, file1, file2, sample, fastqc_bin)
	fastqc_files = get_files(path_to_sample)
	
	for files in fastqc_files:
		parse_fastqcFile(files)
		
############
'''******************************************'''
if __name__== "__main__":
	main()








