#usr/bin/env python
'''
This code calls SPADES assembler and plasmidSPADES.
Retrieves plasmids if any and discard them from the main assembly.
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
from Bio import SeqIO
import shutil

## import my modules
pythonDir = os.path.dirname(os.path.realpath(__file__))
sys.path.append(pythonDir)
import functions

## import configuration
configDir = os.path.dirname(os.path.realpath(__file__)) + '/../../config/'
sys.path.append(configDir)
import config

from blast_parser import parse

perlDir = os.path.dirname(os.path.realpath(__file__)) + '/../perl'
sys.path.append(perlDir)
contig_stats_script = perlDir + '/contig_stats.pl'
rename_seqs_script = perlDir + '/rename_FASTA_seqs.pl'

######
def run_SPADES_plasmid_assembly(path, file1, file2, sample, SPADES_bin, threads):

	print ('+ Running plasmid assembly...')
	name = sample + '_plasmid'
	options = '--plasmid '
	message_return = run_SPADES(path, file1, file2, name, SPADES_bin, options, threads)

	if 	message_return == 'FAIL':	
		print ("\n\n***ERROR: plasmidSPADES failed for sample " + sample)	
		exit()

	scaffolds_retrieved = get_files(path + '/' + name)
	if scaffolds_retrieved == 'FAIL':	
		print ('\n\n***ATTENTION: No plasmids assembly...')

	return (scaffolds_retrieved)

######
def run_SPADES_assembly(path, file1, file2, sample, SPADES_bin, threads):

	print ('+ Running main assembly...')
	options = ''
	message_return = run_SPADES(path, file1, file2, sample, SPADES_bin, options, threads)
	if 	message_return == 'FAIL':	
		print ("\n\n***ERROR: SPADES failed for sample " + sample)
		exit()

	scaffolds_retrieved = get_files(path + '/' + sample)
	if scaffolds_retrieved == 'FAIL':	
		print ('\n\n***ERROR: No scaffolds assembly...')
		exit()
	
	return (scaffolds_retrieved)

######
def run_SPADES(path, file1, file2, name, SPADES_bin, options, threads):

	## call system for SPADES sample given
	sample_folder = functions.create_subfolder(name, path)
	logFile = path + '/' + name + '.log'
	
	if os.path.isfile(logFile):
		return ('OK')
	
	print ("+ Calling spades assembly for sample...", name)	
	cmd_SPADES = '%s %s-t %s -o %s -1 %s -2 %s > %s 2> %s' %(SPADES_bin, options, threads, sample_folder, file1, file2, logFile, logFile)
	print ('[ System Call: ' + cmd_SPADES + ' ]')
	## send command	
	#return (functions.system_call( cmd_SPADES ))
	return ('OK')

######
def get_files(path):
	files = os.listdir(path)
	scaffolds_files = ()
	
	## get files generated with summary information
	for f in files:
		#print (f)
		if f == 'scaffolds.fasta':
			tmp = path + '/' + f 
			#print (tmp)
			
			if os.path.isfile(tmp):
				scaffolds_files = tmp
	if not scaffolds_files:
		return 'FAIL'
	else:
		return scaffolds_files

######
def run_module_SPADES(name, path, file1, file2, threads):

	print ('+ Running modules SPADES...')
	folder = functions.create_subfolder(name, path)
	
	## get configuration
	SPADES_bin = config.EXECUTABLES['spades']
	
	## assembly main 
	path_to_contigs = run_SPADES_assembly(folder, file1, file2, sample, SPADES_bin, threads)

	## assembly plasmids
	path_to_plasmids = run_SPADES_plasmid_assembly(folder, file1, file2, sample, SPADES_bin, threads)
	
	## discard plasmids from main
	tmp_contigs, tmp_plasmids = discardPlasmids(path_to_contigs, path_to_plasmids, folder, sample)
	
	## rename fasta sequences
	new_contigs = tmp_contigs.split("fna.tmp") + '.fna'	
	rename_contigs(tmp_contigs, "scaffolds_chr", new_contigs)
	
	new_plasmids=""
	if os.path.isfile(tmp_plasmids):
		new_plasmids = tmp_contigs.split("fna.tmp") + '.fna'	
		rename_contigs(tmp_plasmids, "scaffolds_plasmids", new_plasmids)
	
	## contig stats
	contig_stats(new_contigs, new_plasmids)	

######
def discardPlasmids(contigs, plasmids, path, sample):
	
	## check if any plasmids
	if (plasmids == 'FAIL'):
		#print ('+ No plasmids assembled.')
		#print ('+ No need to discard any plasmids from the main assembly')
		
		contig_out_file = os.path.dirname(path) + '/' + sample + '/' + sample + '_chromosome.fna.tmp'
		shutil.copy(contigs, contig_out_file)
		return (contig_out_file, plasmids)
	
	## discard 
	print ('+ Check if any plasmids are also reported in main assembly...')

	folder = functions.create_subfolder('blast_search', path)
	
	## generate blastdb for genome
	makeblastDBexe = config.EXECUTABLES['makeblastdb']
	DBname = folder + '/mainAssembly'
	cmd_makeblast = "%s -in %s -input_type fasta -dbtype %s -out %s" %(makeblastDBexe, contigs, 'nucl', DBname)
	print ('[ System Call: ' + cmd_makeblast + ' ]')
	code = functions.system_call(cmd_makeblast)

	if (code == 'FAIL'):
		print ('****ERROR: Some error happened during the makeblastDB command')
		print (cmd_makeblast)
		exit()
	
	## blastn plasmids vs contigs
	blastnexe = config.EXECUTABLES['blastn']
	outFile = folder + '/blastn_output.txt'
	cmd_blastn = "%s -db %s -query %s -out %s -evalue 1e-20 -outfmt \'6 std qlen slen\'" %(blastnexe, DBname, plasmids, outFile )
	print ('[ System Call: ' + cmd_blastn + ' ]')
	codeBlastn = functions.system_call(cmd_blastn)
	
	if (codeBlastn == 'FAIL'):
		print ('****ERROR: Some error happened during the blastn command')
		print (cmd_blastn)
		exit()
	
	########################
	## parseBlast results
	########################
	
	## thresholds
	eval_thresh_float = float(1e-20)
	aln_thresh_given = 90
	min_length = 1000
	
	outFile_parsed = folder + '/blastn_output_parsed.txt'
	output_file = open(outFile_parsed, 'w')	
	sequences2discard = []

	print ('+ Parsing BLAST results generated...\n')
	## get results
	fh = open(outFile)
	for blast_record in parse(fh, eval_thresh=eval_thresh_float, aln_thresh=aln_thresh_given, length_thresh=min_length):
		for hit in blast_record.hits:
			for hsp in hit:
				print('****Alignment****')
				output_file.write('****Alignment****')
				output_file.write('\n')
				
				print('query id: {}'.format(blast_record.qid))
				output_file.write('query id: {}'.format(blast_record.qid))
				output_file.write('\n')
				
				print('sequence: ', hsp.sid)
				sequences2discard.append(hsp.sid)
				output_file.write('sequence: %s' %hsp.sid)
				output_file.write('\n')
				
				print('e value:', hsp.evalue)
				output_file.write('e value: %s' %hsp.evalue)
				output_file.write('\n')
				
				print('aln:', hsp.length)
				output_file.write('aln: %s' %hsp.length)
				output_file.write('\n')
				
				print('qlen:', hsp.qlen, ' [>', min_length, ' bp]')
				output_file.write('qlen: %s [>%s]' %(hsp.qlen, min_length))
				output_file.write('\n')
				
				aln = (int(hsp.qlen)/int(hsp.slen))*100
				print ('aln/slen:', aln, ' [>', aln_thresh_given ,'%]')
				output_file.write('aln/slen: %s [> %s]' %(aln, aln_thresh_given))
				output_file.write('\n\n')

	fh.close()
	output_file.close()
	
	items = len(sequences2discard)
	print ('There are %s sequences to discard from main assembly identified as plasmids' %items)

	## print filtered contigs
	contig_out_file = os.path.dirname(path) + '/' + sample + '/' + sample + '_chromosome.fna.tmp'
	plasmid_out_file = os.path.dirname(path) + '/' + sample + '/' + sample + '_plasmid.fna.tmp'
		
	contig_out_file_handle = open(contig_out_file, 'w')
	plasmid_out_file_handle = open(plasmid_out_file, 'w')
	
	contig_items = SeqIO.parse(contigs, 'fasta')
	for seq in contig_items:
		if seq.id in sequences2discard:
			plasmid_out_file_handle.write(seq.format("fasta"))
			plasmid_out_file_handle.write('\n')
		else:
			contig_out_file_handle.write(seq.format("fasta"))
			contig_out_file_handle.write('\n')

	contig_out_file_handle.close()
	plasmid_out_file_handle.close()	
	return (contig_out_file, plasmid_out_file)

######
def contig_stats(new_contigs, new_plasmids):
	
	## generate contig statistics
	print ('+ Get assembly statistics:...\n')
	print (' + Main assembly:')
	contig_out = new_contigs + '_stats.txt'
	cmd_stats_chromosome = 'perl %s %s > %s' %(contig_stats_script, new_contigs, contig_out)
	code_chr = functions.system_call(cmd_stats_chromosome)
	print (cmd_stats_chromosome)
	
	## dump in screen
	contig_out_file = open(contig_out, 'r')
	contig_out_file_read = contig_out_file.read()
	contig_out_file.close()
	print (contig_out_file_read)
	print ('')	
	if (new_plasmids == 'FAIL'):
		print ('+ No plasmids identified...\n')
	else:
		print ('+ Plasmids assembly')
		plasmid_out = new_plasmids + '_stats.txt'
		cmd_stats_plasmids = 'perl %s %s > %s' %(contig_stats_script, new_plasmids, plasmid_out)
		code_chr = functions.system_call(cmd_stats_plasmids)
		print (cmd_stats_plasmids)

		## dump in screen
		plasmid_out_file = open(plasmid_out, 'r')
		plasmid_file_read = plasmid_out_file.read()
		plasmid_out_file.close()
		print(plasmid_file_read)	
	
######
def rename_contigs(fasta_file, name, new_fasta):
	## perl tools/perl/rename_FASTA_seqs.pl fasta_file name_file name2add ADD|REPLACE|BEGIN|ADD_BEGIN
	perl_call = "%s %s %s %s REPLACE" %(rename_seqs_script, fasta_file, new_fasta, name)
	return (functions.system_call(perl_call))
	
######
def	help_options():
	print ("\nUSAGE: python %s file1 file2 name SPADES_bin threads path\n"  %os.path.realpath(__file__))

######
def main():
	## this code runs when call as a single script

  	## control if options provided or help
	if len(sys.argv) > 1:
		print ("")
	else:
		help_options()
		exit()    	

	file1 = os.path.abspath(argv[1])
	file2 = os.path.abspath(argv[2])
	sample = argv[3]
	SPADES_bin = argv[4]
	threads = int(argv[5])
	path = 	argv[6]

	folder = functions.create_subfolder(sample, path)

	## assembly main 
	path_to_contigs = run_SPADES_assembly(folder, file1, file2, sample, SPADES_bin, threads)
	
	## assembly plasmids
	path_to_plasmids = run_SPADES_plasmid_assembly(folder, file1, file2, sample, SPADES_bin, threads)

	## discard plasmids from main
	tmp_contigs, tmp_plasmids = discardPlasmids(path_to_contigs, path_to_plasmids, folder, sample)
	
	## rename fasta sequences
	new_contigs_list = tmp_contigs.split(".tmp")
	new_contigs = new_contigs_list[0]
	rename_contigs(tmp_contigs, "scaffolds_chr", new_contigs)
	
	new_plasmids=""
	if os.path.isfile(tmp_plasmids):
		new_plasmids_list = tmp_plasmids.split(".tmp")
		new_plasmids = new_plasmids_list[0]
		rename_contigs(tmp_plasmids, "scaffolds_plasmids", new_plasmids)
	
	## contig stats
	contig_stats(new_contigs, new_plasmids)	
		
######

'''******************************************'''
if __name__== "__main__":
	main()


