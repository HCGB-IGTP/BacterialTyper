#!/usr/bin/env python3
##########################################################
## Jose F. Sanchez                                      ##
## Copyright (C) 2019-2020 Lauro Sumoy Lab, IGTP, Spain ##
##########################################################
"""
Retrieves promoter sequences for gene ids from profile analysis generated
"""
## useful imports
from Bio import SeqIO
import sys

#######################################
def get_promoter(file, geneOfInterest, basePairs, sampleName, option, debug=False):
    ## get promoter either from Genbank or fasta file
    if (option == "gbk"):
        get_promoter_gbk(file, geneOfInterest,  basePairs, sampleName, debug)
    elif(option == "fasta"):
        get_promoter_fasta(file, geneOfInterest,  basePairs, sampleName, debug)

#######################################
def get_promoter_fasta(fasta_file, geneOfInterest,  basePairs, sampleName, debug=False):
    print()

#######################################
def get_promoter_gbk(gbf_file, geneOfInterest,  basePairs, sampleName, debug=False):
    """ Parse GenBank file and retrieve the amount of base pairs desired. 
    """
    
    fastaDict={}
    for rec in SeqIO.parse(gbf_file, "genbank"):
    	ID = rec.id
    	SEQ = rec.seq
    
    	## loop through features
    	for feature in rec.features:
    		if feature.type=="gene":
    			qualif = feature.qualifiers
    			for keys, values in qualif.items():
    				#print (keys)
    				#print (values)
    				if values[0]==geneOfInterest:
    					#print (feature)
    					#print (ID)
    					#print (feature.strand)
    					
    					if int(feature.strand) > 0:
    						#print ("Start promoter: " + str(feature.location.nofuzzy_start-int(basePairs)))
    						#print ("Start gene: " + str(feature.location.nofuzzy_start))
    						#print ("End gene: " + str(feature.location.nofuzzy_end))
    						promoter_seq = SEQ[feature.location.nofuzzy_start-int(basePairs):feature.location.nofuzzy_start]
    						#, feature.location.nofuzzy_end           
    					else:
    						#print ("Start promoter: " + str(feature.location.nofuzzy_end+int(basePairs)))
    						#print ("Start gene: " + str(feature.location.nofuzzy_end))
    						#print ("End gene: " + str(feature.location.nofuzzy_start))
    						promoter_seq = SEQ[feature.location.nofuzzy_end : feature.location.nofuzzy_end +int(basePairs)].reverse_complement()
    
    					## print seq
    					id= sampleName + " promoter_" + basePairs + "_" + geneOfInterest
    					fastaDict[id] =promoter_seq 
    
    return(fastaDict)

#######################################
def help_options():
    print ("\nUSAGE: python %s genbank_file gene_id base_pairs sampleName...\n"  %os.path.realpath(__file__))

#######################################
def help_promoter_genes():
    print (colored("\n\n***** TODO: Generate this help message *****\n\n", 'red'))

#######################################
def main():

    ## control if options provided or help
    if len(sys.argv) > 1:
        print ("")
    else:
        help_options()
        exit()
        
    ## arguments
    gbf_file = sys.argv[1]
    geneOfInterest = sys.argv[2]
    basePairs = sys.argv[3]
    sampleName = sys.argv[4]

    ## Debug mode ON
    fastaDict = get_promoter(gbf_file, geneOfInterest,  basePairs, sampleName, True)
    
    print(fastaDict) ## print to file using loop

'''******************************************'''
if __name__== "__main__":
    main()

