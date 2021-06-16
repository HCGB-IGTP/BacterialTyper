#!/usr/bin/env python3
#################################################################
## Jose F. Sanchez                                             ##
## Copyright (C) 2019-2020 Lauro Sumoy Lab, IGTP, Spain        ##
#################################################################
"""
Parses and prepares for phylogenetic analysis
"""
## useful imports
import time
import io
import os
import re
import sys
from sys import argv
from io import open
from termcolor import colored
import pandas as pd
from ete3 import Tree
from Bio import AlignIO
import numpy
import scipy.cluster.hierarchy as sch
from scipy.sparse.csgraph import minimum_spanning_tree
import networkx
import pylab

## import my modules
from BacterialTyper.config import set_config
from BacterialTyper.scripts import min_hash_caller
import HCGB.functions.system_call_functions as HCGB_sys

##################################
def get_snp_distance(aln_file, mode, countGaps, output, Debug):
    """
    
    """

    (D, labeltext) = snp_distance(aln_file, mode, countGaps, Debug)

    ### Write output
    labeltext_string = [str(x) for x in labeltext] ## generat strings, avoid if integers alone
    labeloutname = output + '.labels.txt'
    with open(labeloutname, 'w') as fp:
        fp.write("\n".join(labeltext_string))
        
    ## save matrix in csv format file
    numpy.set_printoptions(suppress=True)
    numpy.savetxt(output + ".csv", D, delimiter=",", fmt='%5.4g', header= ",".join(labeltext_string) )

    ## create minimum spanning tree
        ## https://docs.scipy.org/doc/scipy/reference/generated/scipy.cluster.hierarchy.linkage.html
    Y = sch.linkage(D, method='single', optimal_ordering=True)
    min_hash_caller.get_Newick_tree(Y, D, labeltext, output)
        
        ##D_network = networkx.convert_matrix.from_numpy_matrix(D)
        ##MST = networkx.minimum_spanning_tree(D_network)
        ##print (MST)
        ##fig3 = pylab.figure(figsize=(15, 10))
        ##networkx.draw_networkx(MST, with_labels=True, node_size = 45)
        ##fig3.savefig(output + '.pdf')
        
    return()

##################################
def snp_distance(aln_file, mode, countGaps, Debug):
    """
    Calculate SNP difference between records in alignment
    
    :param aln_file: Alignment file.
    :param mode: Format of aln_file: [nexus,phylip,clustalw,fasta]
    :param countGaps: True/false for counting gaps as differences.
    :param Debug: True/false for debugging messages.
    
    :type aln_file: string
    :type mode: string 
    :type countGaps: bool
    :type Debug: bool
    
    :returns: Numpy matrix with SNP distance matrix and list of labels.
    
    """
    
    ## read alingment
    alignment = AlignIO.read(aln_file, mode)
    
    # build the distance matrix
    D = numpy.zeros([len(alignment), len(alignment)])

    labeltext = []
    
    for i, record in enumerate(alignment):
        for j, record2 in enumerate(alignment):
            snps_diff = 0
            if (record.id != record2.id):
                snps_diff = distance(record.seq, record2.seq, gaps=countGaps)
        
            if Debug:
                print (colored("** DEBUG: record.id, record2.id, snps_diff", 'yellow'))
                print ("%s vs %s : %i" %(record.id, record2.id, snps_diff))
            D[i][j] = snps_diff
            D[j][i] = snps_diff
        
        labeltext.append(record.id)
        
    ## Debug messages
    if Debug:
        print (colored("\n*** DEBUG: SNP distance matrix *****\n", 'red'))
        print ('D:')
        print (D)
        print ('labeltext:')
        print (labeltext)
        print ('Min similarity in matrix: {:.3f}', numpy.min(D)) ## use this to color accordingly

        
    return (D, labeltext)

##################################
def distance(seq1, seq2, gaps=False):
    """
    Get distance between strings
    
    :param seq1: DNA sequence string 1
    :param seq2: DNA sequence string 2
    :param gaps: True/false for counting gaps as differences
    
    :type seq1: string
    :type seq1: string
    :type gaps: bool
    
    :returns: Integer of snps distance
    """
    
    ## check length
    len1 = len(seq1)
    len2 = len(seq2)
    if (len1 != len2):
        print (colored("** ERROR: Length of DNA strings are not the same...", 'red'))
        exit()
    
    ## compare
    snps = 0
    for pos in range (0, min(len1, len2)) :
        if seq1[pos] != seq2[pos]:
            if (seq1[pos] == '-' or seq2[pos]  == '-'):
                if (gaps):
                    snps += 1
            else:
                snps += 1
               
    return(snps)

##################################
def ml_tree(folder, name, threads, output, Debug):
    """
    Create Maximum Likelihood tree reconstruction 
    
    We use IQ-Tree for the versatility and the ability to automatically set parameters. 
    
    :param folder: Snippy-core folder containing results.
    :param name: Name of the analysis.
    :param Debug: True/false for debugging messages
    
    :type folder: string 
    :type name: string
    :type Debug: bool 
    """
    iqtree_exe = set_config.get_exe('iqtree', Debug) 
    bootstrap_number = '1000'
    aln_file = os.path.join(folder, name + '.aln')
    output_log = os.path.join(output, 'iqtree.error.log')
    output_files = os.path.join(output, 'iqtree_' + name)
    
    iqtree_cmd = '%s -s %s -redo --threads-max %s --prefix %s -B %s 2> %s' %(iqtree_exe, aln_file, 
                                                                      threads, output_files, 
                                                                      bootstrap_number, output_log)
    code = HCGB_sys.system_call(iqtree_cmd)
    
    if code == 'OK':
        return ()
    else:
        print ("Some error occurred...")
        return()

    ## raxml is available as conda package
    ## https://anaconda.org/bioconda/raxml
    ## other possibility is FastTree   
    
##################################
def help_options():
    print ("\nUSAGE: python %s file_aln format[nexus,phylip,clustalw,fasta]\n"  %os.path.realpath(__file__))

##################################
def main():
    ## this code runs when call as a single script

    ## control if options provided or help
    if len(sys.argv) > 1:
        print ("")
    else:
        help_options()
        exit()        

    get_snp_distance(argv[1], argv[2], countGaps=True, output="test_phylo", Debug=True)
        
######
if __name__== "__main__":
    main()
    
    
    
    
    
    
    
    
    
    
    
    
    
