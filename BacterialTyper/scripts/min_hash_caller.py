#!/usr/bin/env python3
##########################################################
## Jose F. Sanchez										##
## Copyright (C) 2019 Lauro Sumoy Lab, IGTP, Spain		##
##########################################################
'''
Calls sourmash_ and generates signatures for each genome, clusterizes and generates plots. 

.. include:: ../../links.inc	 
'''
#################################################
## code taken and adapted from: 
##	https://sourmash.readthedocs.io/en/latest/api-example.html
## 	https://github.com/dib-lab/sourmash/blob/master/sourmash/commands.py
#################################################

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
from Bio import Phylo
from io import StringIO

## import modules sourmash
import csv
import shutil
import sourmash
import screed
import numpy
from sourmash import SourmashSignature, save_signatures, load_one_signature
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import pylab
import scipy.cluster.hierarchy as sch

## import my modules
import HCGB
import HCGB.functions.main_functions as HCGB_main
import HCGB.functions.files_functions as HCGB_files

from BacterialTyper.config import set_config
from BacterialTyper.scripts import database_generator

##################################################
def helpMash():
	print (colored("\n\n***** TODO: Generate this help message *****\n\n", 'red'))

##################################################		
def sketch_database(dict_files, folder, Debug, ksize_n, num_sketch):	
	"""Sketch sequence files
	
	This function generates a sourmash index, also called sketch, of the sequences 
	provided in the folder specified.
	
	For speed reasons, we set force=True in add_sequence step to skip over k-mers containing 
	characters other than ACTG, rather than raising an exception.

	:param dict_files: keys are the names of the files and values are the path to the fasta file
	:param folder:
	:param Debug: True/False to print developing messages.
	:param ksize_n: Kmer size value.
	:param num_sketch: Number of sketches to include in the hash signature. 
	
	:type dict_files: Dictionary
	:type folder: string 
	:type Debug: bool
	:type ksize_n: integer
	:type num_sketch: integet
	
	:returns: List of SourmashSignature signatures (siglist) and absolute path files generated (siglist_file).  
	
	
	.. attention:: The code to implement this API function was taken and adapted from: 
	 
		- https://sourmash.readthedocs.io/en/latest/api-example.html
	
		- https://github.com/dib-lab/sourmash/blob/master/sourmash/commands.py
		
	
	.. seealso:: This function depends on sourmash python module (https://sourmash.readthedocs.io/en/latest/). Some functions employed are:
	
		- :func:`sourmash.MinHash`
		
		- :func:`sourmash.SourmashSignature`
		
		- :func:`sourmash.MinHash.add_sequence`
		
		
	.. include:: ../../links.inc	 
	
	"""
	### Default: set as option
	## num_sketch=5000
	## ksize_n=31
	
	minhashes = {}
	for name,g in dict_files.items():
		print ('\t+ Skecthing sample: ', name)
		E = sourmash.MinHash(n=num_sketch, ksize=ksize_n)	## generate hash according to number of sketches and kmer size
		for record in screed.open(g):
			E.add_sequence(record.sequence, True)
		## in add_sequence and for speed reasons, we set force=True to skip over k-mers containing characters other than ACTG, rather than raising an exception.
		minhashes[name]= E
		
	## Debug messages
	if Debug:
		print (colored("\n*** DEBUG: minhashes *****\n", 'red'))
		print (type(minhashes))	
		print (minhashes)

	siglist = []
	siglist_file = []

	### save as signature
	HCGB_files.create_folder(folder)
	for names,hashes in minhashes.items():
		sig1 = SourmashSignature(hashes, name=names)
		outfile_name = folder + '/' + str(names) + '.sig'
		with open(outfile_name, 'wt') as fp:
			save_signatures([sig1], fp)

		siglist_file.append(outfile_name)
		siglist.append(sig1)
	
	return(siglist_file, siglist)		
	
##################################################		
def read_signature(sigfile, ksize_n):
	"""Loads signatures stored in files
	
	Signature generated are stored as MinHashes in JSON format file. Sourmash implements a function to load signatures.
	
	This function here is just a shortcut for the :func:`sourmash.load_one_signature` function in sourmash_.
	
	:param sigfile: Hash signature stored in JSON format
	:param ksize_n: Kmer size.
	
	:type sigfile: string
	:type ksize_n: integer
	
	
	.. seealso:: This function depends on sourmash python module (https://sourmash.readthedocs.io/en/latest/). Some functions employed are:
	
		- :func:`sourmash.load_one_signature`
		
	
	.. include:: ../../links.inc	 
	"""
	
	## code taken and adapted from: https://sourmash.readthedocs.io/en/latest/api-example.html
	#print ('\t+ Loading signature for comparison...')
	sig = load_one_signature(sigfile, ksize=ksize_n, select_moltype='DNA', ignore_md5sum=False)
	return (sig)

##################################################
def compare(siglist, output, Debug):
	#################################################
	## code taken and adapted from: 
	##	https://sourmash.readthedocs.io/en/latest/api-example.html
	## 	https://github.com/dib-lab/sourmash/blob/master/sourmash/commands.py
	#################################################
   
	# build the distance matrix
	D = numpy.zeros([len(siglist), len(siglist)])
	numpy.set_printoptions(precision=3, suppress=True)

	# do all-by-all calculation
	labeltext = []
	for i, E in enumerate(siglist):
		for j, E2 in enumerate(siglist):
			if i < j:
				continue
			similarity = E.similarity(E2, False) ## False -> ignore abundances
			D[i][j] = similarity
			D[j][i] = similarity

		#labeltext.append(E.name())
		labeltext.append(E.name)

   	## Debug messages
	if Debug:
		print (colored("\n*** DEBUG: compare minhashes *****\n", 'red'))
		print (type(D))
		print ('D:')
		print (D)
		print ("Labeltext:")
		print (labeltext)
		print ('Min similarity in matrix: {:.3f}', numpy.min(D)) ## use this to color accordingly
		
	### Write output
	labeltext_string = [str(x) for x in labeltext] ## generat strings, avoid if integers alone
	labeloutname = output + '.labels.txt'
	with open(labeloutname, 'w') as fp:
		fp.write("\n".join(labeltext_string))
		
	## save matrix in csv format file
	numpy.savetxt(output + ".matrix.csv", D, delimiter=",", header= ",".join(labeltext_string) )

	## Debug messages
	if Debug:
		print('saving labels:', labeltext_string)
		print('saving labels to:', labeloutname)
		print('saving distance matrix (binary file) to:', output)
		print('saving distance matrix (csv file) to:', output)

	## return numpy array
	return (D, labeltext_string)

##################################################
def generateNewick(node, newick, parentdist, leaf_names):
	"""Generates Newick tree
	
	This functions generates a newick tree given the a tree object (output from the scipy cluster hierarchy to_tree), the labels and the corresponding matrix distance.
	
	.. note:: scipy.cluster.hierarchy.to_tree
	
		Converts a linkage matrix into an easy-to-use tree object. 
	
	:param node: Tree object (scipy.cluster.hierarchy.to_tree)
	:param newick: Newick tree format
	:param parendist: node.dist or distance matrix
	:param leaf_names: list of labels
	 
	:type newick: string
	:type parendist: matrix
	:type leaf_names: list
	 
	As generateNewick tree is called recursively to create the tree with the parents and children nodes, the first call the argument newick can be empty ("").
	
	.. attention:: The original code was taken and adapted from: 
	
		- https://github.com/scipy/scipy/issues/8274
	
		- https://stackoverflow.com/questions/28222179/save-dendrogram-to-newick-format/31878514#31878514
	
	
	.. include:: ../../links.inc	 
	"""
	
	if node.is_leaf():
		return "%s:%.2f%s" % (leaf_names[node.id], parentdist - node.dist, newick)
	else:
		if len(newick) > 0:
			newick = "):%.2f%s" % (parentdist - node.dist, newick)
		else:
			newick = ");"
		newick = generateNewick(node.get_left(), newick, node.dist, leaf_names)
		newick = generateNewick(node.get_right(), ",%s" % (newick), node.dist, leaf_names)
		newick = "(%s" % (newick)
		return newick

##################################################
def plot(D, labeltext, filename, pdf, colorLabel):
	#################################################
	## code taken and adapted from: 
	##	https://sourmash.readthedocs.io/en/latest/api-example.html
	## 	https://github.com/dib-lab/sourmash/blob/master/sourmash/commands.py
	#################################################

	matrix_out = filename + '.matrix'

	### decide on PDF/PNG output
	if pdf:
		matrix_out += '.pdf'
	
	else:
		matrix_out += '.png'
	
	###########################
	### make the histogram
	###########################
	#print ('+ Saving histogram of matrix values: ', hist_out)
	#fig = pylab.figure(figsize=(8,5))
	#pylab.hist(numpy.array(D.flat), bins=100)
	#fig.savefig(hist_out)

	#######################################
	### make the dendrogram: do clustering
	## https://docs.scipy.org/doc/scipy/reference/cluster.hierarchy.html
	#######################################
	#fig2 = pylab.figure(figsize=(11,8))
	#ax1 = fig2.add_axes([0.1, 0.1, 0.7, 0.8])
	#ax1.set_xticks([])
	#ax1.set_yticks([])
	#Y = sch.linkage(D, method='single')
	#Z1 = sch.dendrogram(Y, orientation='right', labels=labeltext)
	#fig2.savefig(dendrogram_out)
	#print ('+ Wrote dendrogram to:', dendrogram_out)

	#######################################
	### make the dendrogram+matrix:
		#######################################

	### Original code
	## fig = sourmash_fig.plot_composite_matrix(D, labeltext, show_labels=args.labels,
	##			 show_indices=args.indices, vmin=args.vmin, vmax=args.vmax, force=args.force)
	## I get code from the source code for this function and use it here.
	## Get to generate a slightly different image representation

	Y = sch.linkage(D, method='single') ## https://docs.scipy.org/doc/scipy/reference/generated/scipy.cluster.hierarchy.linkage.html
	fig3 = pylab.figure(figsize=(15, 10))
	ax1 = fig3.add_axes([0.09, 0.1, 0.2, 0.6])

	# plot dendrogram
	## https://python-graph-gallery.com/400-basic-dendrogram/
	Z1_2 = sch.dendrogram(Y, orientation='left', labels=labeltext)

	## set label colors if desired
	if colorLabel:
		sch.dendrogram(Y, orientation='left', labels=labeltext)
		ax = plt.gca()
		xlbls = ax.get_ymajorticklabels()
		for lbl in xlbls:
			lbl.set_color(colorLabel[ lbl.get_text() ])

	ax1.set_xticks([])
	xstart = 0.45
	width = 0.45
	scale_xstart = xstart + width + 0.01

	# plot matrix
	axmatrix = fig3.add_axes([xstart, 0.1, width, 0.6])

	# (this reorders D by the clustering in Z1_2)
	idx1 = Z1_2['leaves']
	D = D[idx1, :]
	D = D[:, idx1]

 	## use minimun distance matrix D to color accordingly
	min_D = numpy.min(D) ## 
	cmap_palette = ""
	if min_D > 0.8:
		cmap_palette = pylab.cm.YlGnBu
	else:
		cmap_palette = pylab.cm.YlGnBu

	# show matrix
	im = axmatrix.matshow(D, aspect='auto', origin='lower', cmap= cmap_palette, vmin=1, vmax=0)
	axmatrix.set_xticks([])
	axmatrix.set_yticks([])
	
	## https://matplotlib.org/3.1.0/api/_as_gen/matplotlib.pyplot.matshow.html
	## https://matplotlib.org/3.1.0/api/_as_gen/matplotlib.axes.Axes.imshow.html#matplotlib.axes.Axes.imshow
	## cmap=pylab.cm.RdBu
	## https://scipy-cookbook.readthedocs.io/items/Matplotlib_Show_colormaps.html

	# Plot colorbar.
	axcolor = fig3.add_axes([scale_xstart, 0.1, 0.02, 0.6])
	pylab.colorbar(im, cax=axcolor)	
	fig3.savefig(matrix_out)
	print ('+ Wrote matrix to:', matrix_out)
	return(Y)
	
##################################################

def get_Newick_tree(cluster_hierachy, DataMatrix, labeltext, output):
	"""
	
	.. seealso:: This function depends on other BacterialTyper functions called:
	
		- :func:`BacterialTyper.scripts.functions.printList2file`
		
	"""
	
	tree = sch.to_tree(cluster_hierachy, False)
	Newick_tree = generateNewick(tree, "", DataMatrix, labeltext)
	
	Newick_tree_file = output + '.nwk'
	HCGB_main.printList2file(Newick_tree_file, [Newick_tree])
	
	handle = StringIO(Newick_tree)
	treePhylo = Phylo.read(handle, "newick")
	
	leaves_tree = []
	for leaf in treePhylo.get_terminals(): 
		leaves_tree.append(leaf.name)
		#print(leaf.name)
	
	## BUG: the list is printed alphabetically ordered
	Newick_tree_leaves =  output + '.leaves.txt'
	HCGB_main.printList2file(Newick_tree_leaves, leaves_tree)

def	help_options():
	print ("\nUSAGE: python %s fasta_1,ID1 fasta_2,ID2 ... fasta_n,IDn \n"  %os.path.realpath(__file__))

##################################################
def main():
	## this code runs when call as a single script

	## control if options provided or help
	if len(sys.argv) > 1:
		print ("")
	else:
		help_options()
		exit()
	
	file_names_dict = {}
	for i in argv:
		if i == argv[0]:
			continue
		
		file_names_dict[ i.split(',')[1] ] = i.split(',')[0]
	
	## to do: implement main function
	folder = "."
	Debug = False
	output = "test"
	pdf = True
	
	## Defaults
	ksize_n = 51
	num_sketch = 5000
	
	### create signatures and compare them
	(siglist_file, siglist) = sketch_database(file_names_dict, folder, Debug, ksize_n, num_sketch)
	(D, labeltext) = compare(siglist, output, Debug)
	
	## create plot
	Y = plot(D, labeltext, output, pdf, "")
	
	## create newick tree file
	get_Newick_tree(Y, D, labeltext, output)
	

##################################################
if __name__== "__main__":
	main()

	
##################################################		
##def sketch_database(list_files, mash_bin, out_file, name, folder):	
##	### -p threads
##	mash_cmd = ""
##	if len(list_files) == 1:
##		print ('\t+ Skecthing sample: ', name)
##		mash_cmd = '%s sketch -o %s %s' %(mash_bin, out_file, list_files[0])
##
##	else:
##		print ('\t+ Skecthing list of samples provided: ', name)
##		## batch
##		out_batch = folder + "/.batch_entries.txt"
##		functions.printList2file(out_batch, list_files)
##		mash_cmd = '%s sketch -l -o %s %s' %(mash_bin, out_file, out_batch)
##	
##	return(functions.system_call(mash_cmd))

# subsample?
#if args.subsample:
#	numpy.random.seed(args.subsample_seed)

#	sample_idx = list(range(len(labeltext)))
#	numpy.random.shuffle(sample_idx)
#	sample_idx = sample_idx[:args.subsample]

#	np_idx = numpy.array(sample_idx)
#	D = D[numpy.ix_(np_idx, np_idx)]
#	labeltext = [ labeltext[idx] for idx in sample_idx ]
