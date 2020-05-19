#!/usr/bin/env python3
#########################################################
## Jose F. Sanchez                                      ##
## Copyright (C) 2019 Lauro Sumoy Lab, IGTP, Spain      ##
##########################################################
"""
Main script for the BacterialTyper project.
"""
import argparse 
import os
import sys
import BacterialTyper.modules
 

## initiate parser
parser = argparse.ArgumentParser(
    prog='BacterialTyper',
    description='BacterialTyper: Bacterial Genotyping using NGS'
    #,epilog="(c) 2019. Jose F. Sanchez, Cristina Prat and Lauro Sumoy."
    
)
subparsers = parser.add_subparsers(title='Available modules', help='', metavar='')

## help options list
help_options = ('--help_Database',
                '--help_format',
				'--help_BUSCO',
				'--help_project',
				'--help_Prokka',
				'--help_ARIBA',
				'--help_PhiSpy',
				'--help_trimm_adapters',
				'--help_multiqc',
				'--help_KMA',
				'--help_MLSTar',
				'--help_MGE_analysis',
				'--help_Mash',
				'--help_input_MGE',
                '--help_Snippy')

## space
##subparser_space = subparsers.add_parser('-----------------------', help='')
##subparser_space = subparsers.add_parser('Pipeline configuration', help='')
##subparser_space = subparsers.add_parser('-----------------------', help='')
##subparser_space = subparsers.add_parser(' ', help='')

#######################
#### Configuration ####
#######################
## submodules
##------------------------------ config ---------------------- ##

## [TODO]
subparser_config = subparsers.add_parser(
    'config',
    help='Configure the pipeline',
    description='Configure dependencies, executables and additional python or perl modules.',
)
subparser_config.add_argument("option", help="Checks if missing any dependencies or modules or tries to install them.", choices=['check','install'])
subparser_config.add_argument("--install_path", help="Path to install missing modules or dependencies. [Default: BacterialTyper config folder]")
subparser_config.add_argument("--IslandPath", action="store_true", help="Check for additional perl and software packages required for IslandPath.")
subparser_config.add_argument("--debug", action="store_true", help="Show additional message for debugging purposes.")
subparser_config.set_defaults(func=BacterialTyper.modules.config.run)

## add fake module blank to add space
subparser_space = subparsers.add_parser(' ', help='')
##-------------------------------------------------------------##

###################
#### Databases ####
###################

##---------------------  database -------------------- ##
subparser_database = subparsers.add_parser(
    'database',
    help='Initiates/updates a database for later use.',
    description='This module initiates or updates a database with information from multiple sources: NCBI, ARIBA databases, KMA indexes, BUSCO, PubMLST...',
)

initdb_input_options_group = subparser_database.add_argument_group("Input")
initdb_input_options_group.add_argument("--path", help="Folder path to generate database.", required= not any(elem in help_options for elem in sys.argv) )
subparser_database.add_argument("--threads", type=int, help="Number of CPUs to use [Default: 2].", default=2) ## for ARIBA prepareref

## NCBI
initdb_NCBIoptions_group = subparser_database.add_argument_group("NCBI Data")
initdb_NCBIoptions_group.add_argument("--ID_file", help="CSV file containing several columns per row according to the header provided: ##genus,species,name,NCBI_assembly_ID")
initdb_NCBIoptions_group.add_argument("--descendant", help="Get related indexed genomes for a given NCBI taxonomy id", type=int)

## user_data
initdb_user_data_options_group = subparser_database.add_argument_group("Project data")
initdb_user_data_options_group.add_argument("--project_folder", help="Folder containing previously identified samples under mode --project. For more information check --help_Database or --help_project option")
initdb_user_data_options_group.add_argument("--single_end", action="store_true", help="Single end files [Default OFF]. Default mode is paired-end.")
initdb_user_data_options_group.add_argument("--batch", action="store_true", help="Provide this option if input is a file containing multiple paths instead a path.")
initdb_user_data_options_group.add_argument("--in_sample", help="File containing a list of samples to include (one per line) from input folder(s) [Default OFF].")
initdb_user_data_options_group.add_argument("--ex_sample", help="File containing a list of samples to exclude (one per line) from input folder(s) [Default OFF].")

## KMA
initdb_KMAoptions_group = subparser_database.add_argument_group("KMA Databases")
initdb_KMAoptions_group.add_argument("--kma_db", dest='kma_dbs', nargs='*', help="kma database(s) to download. Provide several input if desired. [Default: bacteria & plasmids]", choices=['bacteria', 'archaea', 'protozoa', 'fungi', 'plasmids', 'typestrains'])
initdb_KMAoptions_group.add_argument("--no_def_kma", action="store_true", help="Only applicable if kma_db is ON. It discards default databases (bacteria & plasmids) from kma_db. Only user selected databases indexed. [Default: OFF]")
##initdb_KMAoptions_group.add_argument("--index_kma", action="store_true", help="Index the genomes downloaded for later usage during species identification. Not compatible with --kma_db options: available databases with broader taxonomic ranges. [Default: OFF]")

## ariba
initdb_ARIBAoptions_group = subparser_database.add_argument_group("ARIBA Databases")
initdb_ARIBAoptions_group.add_argument("--ARIBA_db", dest='ariba_dbs', nargs='*', help="ARIBA database(s) to download. Provide several input if desired. Not compatible with --no_ARIBA. [Default: CARD & VFDB]", choices=['CARD', 'MEGARes', 'VFDB', 'ResFinder', 'PlasmidFinder', 'VirulenceFinder', 'srst2', 'ARG-ANNOT'])
initdb_ARIBAoptions_group.add_argument("--no_def_ARIBA", action="store_true", help="It discards default databases (CARD & VFDB) from ARIBA_db. Only user selected databases indexed. [Default: OFF]")
initdb_ARIBAoptions_group.add_argument("--user_ARIBA_db_fasta", dest='ariba_users_fasta', nargs='*', help="User selected database: fasta sequences. [Default: OFF]")
initdb_ARIBAoptions_group.add_argument("--user_ARIBA_db_metadata", dest='ariba_users_meta', nargs='*', help="User selected database: metadata. [Default: OFF]") 
initdb_ARIBAoptions_group.add_argument("--no_ARIBA", action="store_true", help="It prevents downloading ARIBA databases. [Default: OFF]")

## BUSCO 
initdb_BUSCOoptions_group = subparser_database.add_argument_group("BUSCO Datasets")
initdb_BUSCOoptions_group.add_argument("--BUSCO_dataset", dest='BUSCO_dbs', nargs='*', help="BUSCO dataset to use according to your sample(s) expected taxonomic range. Provide several input if desired. For more information provide --help_BUSCO option")

## help messages
info_group_database = subparser_database.add_argument_group("Additional information")
info_group_database.add_argument("--debug", action="store_true", help="Show additional message for debugging purposes.")
info_group_database.add_argument("--help_Database", action="store_true", help="Print further information for this database module")
info_group_database.add_argument("--help_ARIBA", action="store_true", help="Print further information for ARIBA databases")
info_group_database.add_argument("--help_BUSCO", action="store_true", help="Benchmarking Universal Single Copy Orthologs (BUSCO) dataset help information.")
info_group_database.add_argument("--help_KMA", action="store_true", help="Show additional help on KMA software and options.")
info_group_database.add_argument("--help_project", action="store_true", help="Show additional help on the project scheme.")

subparser_database.set_defaults(func=BacterialTyper.modules.database.run_database)
##-------------------------------------------------------------##

## space
subparser_space = subparsers.add_parser(' ', help='')

#########################
#### Prepare samples ####
#########################

##--------------------------- prepareSamples ----------------- ##
subparser_prep = subparsers.add_parser(
    'prep',
    help='Prepares FASTQ files from samples',
    description='This module prepares fastq files from a sequencing run. It could renamed, copy, link or merge them when multiples files have been generated for the same sample e.g different lanes. It concatenates these files according the common identifier and generates a unique file, one per paired-read if necessary',
)

prep_mode_name = subparser_prep.add_argument_group("Mode")
prep_mode = prep_mode_name.add_mutually_exclusive_group(required= not any(elem in help_options for elem in sys.argv))
prep_mode.add_argument("--project", action="store_true", help="It initiates --output_folder folder to contain a project with samples, metadata, configuration etc. [Default]")
prep_mode.add_argument("--detached", action="store_true", help="Isolated mode. No project folder initiated for further steps.")

in_out_group_prep = subparser_prep.add_argument_group("Input/Output")
in_out_group_prep.add_argument("--input", help="Folder containing fastq files. Files could be .fastq/.fq/ or fastq.gz/.fq.gz. All files would be retrieved.", required= not any(elem in help_options for elem in sys.argv))
in_out_group_prep.add_argument("--output_folder", help="Output folder. Name for the project folder.", required= not any(elem in help_options for elem in sys.argv))
in_out_group_prep.add_argument("--single_end", action="store_true", help="Single end files [Default OFF]. Default mode is paired-end.")
in_out_group_prep.add_argument("--batch", action="store_true", help="Provide this option if input is a file containing multiple paths instead a path.")
in_out_group_prep.add_argument("--in_sample", help="File containing a list of samples to include (one per line) from input folder(s) [Default OFF].")
in_out_group_prep.add_argument("--ex_sample", help="File containing a list of samples to exclude (one per line) from input folder(s) [Default OFF].")

options_group_prep = subparser_prep.add_argument_group("Options")
options_group_prep.add_argument("--threads", type=int, help="Number of CPUs to use [Default: 2].", default=2)
options_group_prep.add_argument("--copy", action="store_true", help="Instead of generating symbolic links, copy files into output folder. [Default OFF].")
options_group_prep.add_argument("--merge", action="store_true", help="Merges FASTQ files for the same sample [Default OFF].")
options_group_prep.add_argument("--rename", help="File containing original name and final name for each sample separated by comma. No need to provide a name for each pair if paired-end files. If provided with option '--merge', the merge files would be renamed accordingly.")

info_group_prep = subparser_prep.add_argument_group("Additional information")
info_group_prep.add_argument("--help_format", action="store_true", help="Show additional help on name format for files.")
info_group_prep.add_argument("--debug", action="store_true", help="Show additional message for debugging purposes.")

subparser_prep.set_defaults(func=BacterialTyper.modules.prep.run_prep)
##-------------------------------------------------------------##

##--------------------------- QC ------------------------- ##
subparser_qc = subparsers.add_parser(
    'QC',
    help='Quality check for samples',
    description='This module calls different programs attending the input provided: FASTQC for a quality check of raw reads; BUSCO for assembly and/or annotations provided for sample.',
)
## other options
annotation_options = ('--assembly', '--annotation')
detached_options = ('--detached', '--output_folder')

qc_mode_name = subparser_qc.add_argument_group("Mode")
qc_mode = qc_mode_name.add_mutually_exclusive_group(required= not any(elem in help_options for elem in sys.argv))
qc_mode.add_argument("--project", action="store_true", help="Requires as --input a folder containing a project with samples, metadata, configuration etc. [Default]")
qc_mode.add_argument("--detached", action="store_true", help="Isolated mode. --input is a folder containing samples, contigs or protein sequences. Provide a unique path o several using --batch option")

in_out_group_qc = subparser_qc.add_argument_group("Input/Output")
in_out_group_qc.add_argument("--input", help="Folder containing input. Project or raw reads, assembly or annotation fasta files according to mode option provided.", required= not any(elem in help_options for elem in sys.argv))
in_out_group_qc.add_argument("--output_folder", help="Output folder. Required if '--detached' mode. Under '--project' mode, information will be stored following a designed scheme. See instructions for further details", required = '--detached' in sys.argv)
in_out_group_qc.add_argument("--batch", action="store_true", help="Provide this option if input is a file containing multiple paths instead a path.")
in_out_group_qc.add_argument("--in_sample", help="File containing a list of samples to include (one per line) from input folder(s) [Default OFF].")
in_out_group_qc.add_argument("--ex_sample", help="File containing a list of samples to exclude (one per line) from input folder(s) [Default OFF].")

exclusive_group_qc_name = subparser_qc.add_argument_group("Options")
exclusive_group_qc = exclusive_group_qc_name.add_mutually_exclusive_group(required= not any(elem in help_options for elem in sys.argv))
exclusive_group_qc.add_argument("--raw_reads", action="store_true",  help="Check quality for each sample using FASTQC analysis. Input: reads (fastq/fq). See --help_format for further details.")
exclusive_group_qc.add_argument("--assembly", action="store_true",  help="Check assembly completeness using BUSCO and descriptive statistics. Input: draft assemblies, scaffolds, contigs... See --help_BUSCO for additional details.")
exclusive_group_qc.add_argument("--annotation", action="store_true",  help="Check annotation completenes using BUSCO statistics. Input: protein sequences in fasta format. See --help_BUSCO for additional details.")

options_group_qc = subparser_qc.add_argument_group("Configuration")
options_group_qc.add_argument("--single_end", action="store_true", help="Single end files [Default OFF]. Default mode is paired-end. Only applicable if --raw_reads option.")
options_group_qc.add_argument("--skip_report", action="store_true", help="Do not report statistics using MultiQC report module [Default OFF]")
options_group_qc.add_argument("--threads", type=int, help="Number of CPUs to use [Default: 2].", default=2)

dataset_group_qc = subparser_qc.add_argument_group("Datasets")
dataset_group_qc.add_argument("--database", help="Directory containing databases previously downloaded such as ARIBA, KMA, BUSCO genbank and user_data folders.", required= any(elem in annotation_options for elem in sys.argv))
dataset_group_qc.add_argument("--BUSCO_dataset", dest='BUSCO_dbs', nargs='*', help="BUSCO dataset to use according to your sample(s) expected taxonomic range. Provide several input if desired. For more information provide --help_BUSCO option", required= any(elem in annotation_options for elem in sys.argv))

info_group_qc = subparser_qc.add_argument_group("Additional information")
info_group_qc.add_argument("--help_BUSCO", action="store_true", help="Benchmarking Universal Single Copy Orthologs (BUSCO) dataset help information.")
info_group_qc.add_argument("--debug", action="store_true", help="Show additional message for debugging purposes.")
info_group_qc.add_argument("--help_format", action="store_true", help="Show additional help on name format for files.")
info_group_qc.add_argument("--help_project", action="store_true", help="Show additional help on the project scheme.")
info_group_qc.add_argument("--help_multiqc", action="store_true", help="Show additional help on the multiQC module.")
subparser_qc.set_defaults(func=BacterialTyper.modules.qc.run_QC)
##-------------------------------------------------------------##

##------------------------------ trimm ----------------------- ##
subparser_trimm = subparsers.add_parser(
    'trimm',
    help='Trimms sequencing adapters.',
    description='This module trimms sequencing adapters that could be present in next generation sequencing files',
)
trimm_mode_name = subparser_trimm.add_argument_group("Mode")
trimm_mode = trimm_mode_name.add_mutually_exclusive_group(required= not any(elem in help_options for elem in sys.argv))
trimm_mode.add_argument("--project", action="store_true", help="Requires as --input a folder containing a project with samples, metadata, configuration etc. [Default]")
trimm_mode.add_argument("--detached", action="store_true", help="Isolated mode. --input is a folder containing fastq reads. Provide a unique path o several using --batch option")

in_out_group_trimm = subparser_trimm.add_argument_group("Input/Output")
in_out_group_trimm.add_argument("--input", help="Folder containing a project or reads, according to the mode selected. Files could be .fastq/.fq/ or fastq.gz/.fq.gz. See --help_format for additional details.", required= not any(elem in help_options for elem in sys.argv))
in_out_group_trimm.add_argument("--output_folder", help="Output folder.", required = '--detached' in sys.argv)
in_out_group_trimm.add_argument("--single_end", action="store_true", help="Single end files [Default OFF]. Default mode is paired-end.")
in_out_group_trimm.add_argument("--batch", action="store_true", help="Provide this option if input is a file containing multiple paths instead a path.")
in_out_group_trimm.add_argument("--in_sample", help="File containing a list of samples to include (one per line) from input folder(s) [Default OFF].")
in_out_group_trimm.add_argument("--ex_sample", help="File containing a list of samples to exclude (one per line) from input folder(s) [Default OFF].")

options_group_trimm = subparser_trimm.add_argument_group("Options")
options_group_trimm.add_argument("--skip_report", action="store_true", help="Do not report statistics using MultiQC report module [Default OFF]. See details in --help_multiqc")
options_group_trimm.add_argument("--adapters", help="Adapter sequences to use for the trimming process. See --help_trimm_adapters for further information.")
options_group_trimm.add_argument("--threads", type=int, help="Number of CPUs to use [Default: 2].", default=2)

info_group_trimm = subparser_trimm.add_argument_group("Additional information")
info_group_trimm.add_argument("--debug", action="store_true", help="Show additional message for debugging purposes.")
info_group_trimm.add_argument("--help_format", action="store_true", help="Show additional help on name format for files.")
info_group_trimm.add_argument("--help_trimm_adapters", action="store_true", help="Show additional information on trimm adapters.")
info_group_trimm.add_argument("--help_project", action="store_true", help="Show additional help on the project scheme.")
info_group_trimm.add_argument("--help_multiqc", action="store_true", help="Show additional help on the multiQC module.")

subparser_trimm.set_defaults(func=BacterialTyper.modules.trimm.run)
##-------------------------------------------------------------##

##################
#### Assembly ####
##################
##--------------------------- assemble ----------------------- ##
subparser_assemble = subparsers.add_parser(
    'assemble',
    help='Assembly for each sample.',
    description='This module assembles sequencing reads into contigs using spades. Input must be trimmed sequences.'
)
assembly_mode_name = subparser_assemble.add_argument_group("Mode")
assembly_mode = assembly_mode_name.add_mutually_exclusive_group(required= not any(elem in help_options for elem in sys.argv))
assembly_mode.add_argument("--project", action="store_true", help="Requires as --input a folder containing a project with samples, metadata, configuration etc. [Default]")
assembly_mode.add_argument("--detached", action="store_true", help="Isolated mode. --input is a folder containing trimmed fastq reads. Provide a unique path o several using --batch option")

in_out_group_assembly = subparser_assemble.add_argument_group("Input/Output")
in_out_group_assembly.add_argument("--input", help="Folder containing a project or reads, according to the mode selected. Files could be .fastq/.fq/ or fastq.gz/.fq.gz. See --help_format for additional details.", required= not any(elem in help_options for elem in sys.argv))
in_out_group_assembly.add_argument("--output_folder", help="Output folder.", required = '--detached' in sys.argv)
in_out_group_assembly.add_argument("--single_end", action="store_true", help="Single end files [Default OFF]. Default mode is paired-end.")
in_out_group_assembly.add_argument("--batch", action="store_true", help="Provide this option if input is a file containing multiple paths instead a path.")
in_out_group_assembly.add_argument("--in_sample", help="File containing a list of samples to include (one per line) from input folder(s) [Default OFF].")
in_out_group_assembly.add_argument("--ex_sample", help="File containing a list of samples to exclude (one per line) from input folder(s) [Default OFF].")

dataset_group_assembly = subparser_assemble.add_argument_group("Datasets")
dataset_group_assembly.add_argument("--database", help="Directory containing databases previously downloaded such as ARIBA, KMA, BUSCO genbank and user_data folders.", required=not any(elem in help_options for elem in sys.argv))
dataset_group_assembly.add_argument("--BUSCO_dataset", dest='BUSCO_dbs', nargs='*', help="BUSCO dataset to use according to your sample(s) expected taxonomic range. Provide several input if desired. For more information provide --help_BUSCO option", required= not any(elem in help_options for elem in sys.argv))

options_group_assembly = subparser_assemble.add_argument_group("Options")
options_group_assembly.add_argument("--threads", type=int, help="Number of CPUs to use [Default: 2].", default=2)
dataset_group_assembly.add_argument("--skip_report", action="store_true", help="Do not report statistics using MultiQC report module [Default OFF]. See details in --help_multiqc")

info_group_assemble = subparser_assemble.add_argument_group("Additional information")
info_group_assemble.add_argument("--debug", action="store_true", help="Show additional message for debugging purposes.")
info_group_assemble.add_argument("--help_multiqc", action="store_true", help="Show additional help on the multiQC module.")
info_group_assemble.add_argument("--help_format", action="store_true", help="Show additional help on name format for files.")
info_group_assemble.add_argument("--help_project", action="store_true", help="Show additional help on the project scheme.")
info_group_assemble.add_argument("--help_BUSCO", action="store_true", help="Benchmarking Universal Single Copy Orthologs (BUSCO) dataset help information.")

subparser_assemble.set_defaults(func=BacterialTyper.modules.assemble.run_assembly)
##-------------------------------------------------------------##

####################
#### Annotation ####
####################
##--------------------------- annotate ----------------------- ##
subparser_annotate = subparsers.add_parser(
    'annot',
    help='Annotation for each sample.',
    description='This module annotates contig/scaffold assemblies and generates protein, gff and other annotation information. Input must be fasta/fna assemblies.'
)
annotate_mode_name = subparser_annotate.add_argument_group("Mode")
annotate_mode = annotate_mode_name.add_mutually_exclusive_group(required= not any(elem in help_options for elem in sys.argv))
annotate_mode.add_argument("--project", action="store_true", help="Requires as --input a folder containing a project with samples, metadata, configuration etc. [Default]")
annotate_mode.add_argument("--detached", action="store_true", help="Isolated mode. --input is a folder containing assemblies or scaffolds. Provide a unique path o several using --batch option")

in_out_group_annot = subparser_annotate.add_argument_group("Input/Output")
in_out_group_annot.add_argument("--input", help="Folder containing a project or assemblies. See --help_format for additional details.", required= not any(elem in help_options for elem in sys.argv) )
in_out_group_annot.add_argument("--output_folder", help="Output folder.", required = '--detached' in sys.argv)
##in_out_group_annot.add_argument("--batch", help="Provide a csv file containing the name and the path for each assembly. No header. Provided it in format: name,tag,file. tag = chromosome/plasmid. e.g. sample1,chromosome,/path/to/sample1/assembly/file.fasta\nsample1,plasmid,/path/to/sample1/assembly/file.fasta")
in_out_group_annot.add_argument("--in_sample", help="File containing a list of samples to include (one per line) from input folder(s) [Default OFF].")
in_out_group_annot.add_argument("--ex_sample", help="File containing a list of samples to exclude (one per line) from input folder(s) [Default OFF].")

dataset_group_annot = subparser_annotate.add_argument_group("Datasets")
dataset_group_annot.add_argument("--database", help="Directory containing databases previously downloaded such as ARIBA, KMA, BUSCO genbank and user_data folders.", required= not any(elem in help_options for elem in sys.argv) )
dataset_group_annot.add_argument("--BUSCO_dataset", dest='BUSCO_dbs', nargs='*', help="BUSCO dataset to use according to your sample(s) expected taxonomic range. Provide several input if desired. For more information provide --help_BUSCO option", required= not any(elem in help_options for elem in sys.argv) )

param_group_annot = subparser_annotate.add_argument_group("Parameters")
param_group_annot.add_argument("--skip_report", action="store_true", help="Do not report statistics using MultiQC report module [Default OFF]")
param_group_annot.add_argument("--threads", type=int, help="Number of CPUs to use [Default: 2].", default=2)

prokka_group_annot = subparser_annotate.add_argument_group("Prokka Options")
prokka_group_annot.add_argument("--kingdom", help="Select a kingdom fitting your samples.", choices=['Archaea','Bacteria','Mitochondria','Viruses'], required= not any(elem in help_options for elem in sys.argv) )
prokka_group_annot.add_argument("--genera", help="Select a genera fitting your samples.", choices=['Enterococcus','Escherichia','Staphylococcus', 'Other'], required= not any(elem in help_options for elem in sys.argv) )

info_group_annot = subparser_annotate.add_argument_group("Additional information")
info_group_annot.add_argument("--debug", action="store_true", help="Show additional message for debugging purposes.")
info_group_annot.add_argument("--help_BUSCO", action="store_true", help="Benchmarking Universal Single Copy Orthologs (BUSCO) dataset help information.")
info_group_annot.add_argument("--help_Prokka", action="store_true", help="Prokkaryotic annotation help information.")
info_group_annot.add_argument("--help_multiqc", action="store_true", help="Show additional help on the multiQC module.")
info_group_annot.add_argument("--help_format", action="store_true", help="Show additional help on name format for files.")
info_group_annot.add_argument("--help_project", action="store_true", help="Show additional help on the project scheme.")

subparser_annotate.set_defaults(func=BacterialTyper.modules.annot.run_annotation)
##-------------------------------------------------------------##

## space
subparser_space = subparsers.add_parser(' ', help='')

########################
#### Identification ####
########################

##--------------------------- identification ------------------##
subparser_ident = subparsers.add_parser(
    'ident',
    help='Species identification for each sample.',
    description='This module calls a kmer strategy for a species and strain identification',
)
ident_mode_name = subparser_ident.add_argument_group("Mode")
ident_mode = ident_mode_name.add_mutually_exclusive_group(required= not any(elem in help_options for elem in sys.argv))
ident_mode.add_argument("--project", action="store_true", help="Requires as --input a folder containing a project with samples, metadata, configuration etc. [Default]")
ident_mode.add_argument("--detached", action="store_true", help="Isolated mode. --input is a folder containing assemblies or scaffolds. Provide a unique path o several using --batch option")

initial_group_ident = subparser_ident.add_argument_group("Input/Output")
initial_group_ident.add_argument("--input", help="Folder containing trimmed reads. Files could be .fastq/.fq/ or fastq.gz/.fq.gz", required = not any(elem in help_options for elem in sys.argv))
initial_group_ident.add_argument("--output_folder", help="Output folder.", required = '--detached' in sys.argv)
initial_group_ident.add_argument("--single_end", action="store_true", help="Single end files [Default OFF]. Default mode is paired-end.")
initial_group_ident.add_argument("--batch", action="store_true", help="Provide this option if input is a file containing multiple paths instead a path.")
initial_group_ident.add_argument("--in_sample", help="File containing a list of samples to include (one per line) from input folder(s) [Default OFF].")
initial_group_ident.add_argument("--ex_sample", help="File containing a list of samples to exclude (one per line) from input folder(s) [Default OFF].")

dataset_group_ident = subparser_ident.add_argument_group("Datasets")
dataset_group_ident.add_argument("--database", help="Directory containing databases previously downloaded such as ARIBA, KMA, BUSCO genbank and user_data folders.", required=not any(elem in help_options for elem in sys.argv))
dataset_group_ident.add_argument("--MLST_profile", help="Directory containing MLST profile. By default, a folder named PubMLST is generated under database folder. Provide this option if profile to be used is from another source than PubMLST.org.")

kma_group = subparser_ident.add_argument_group("KMA Databases Configuration")
kma_group.add_argument("--kma_db", dest='kma_dbs', nargs='*', help="kma database(s) to check. Provide several input if desired. If it is not available it would be downloaded and included in default database folder. [Default: bacteria & plasmids]", choices=['bacteria', 'archaea', 'protozoa', 'fungi', 'plasmids', 'typestrains'])
kma_group.add_argument("--kma_external_file", dest='kma_external_files', nargs='*', help="External fasta file to include in the search. It could be indexed or not. [Default OFF].")
kma_group.add_argument("--only_kma_db", action="store_true", help="Provide this option if only user defined kma database(s) provided via --kma_db should be used for species identification. Discard Bacteria and Plasmids as default.")

options_group_name_ident = subparser_ident.add_argument_group("Options")
options_group_name_ident.add_argument("--user_data", action="store_true", help="Use previously identified samples in user_data [Default OFF].")
options_group_name_ident.add_argument("--genbank_data", action="store_true", help="Use reference genomes from NCBI stored in database folder [Default OFF].")

exclusive_group_name_ident2 = subparser_ident.add_argument_group("Exclusive")
exclusive_group2 = exclusive_group_name_ident2.add_mutually_exclusive_group()
exclusive_group2.add_argument("--all_data", action="store_true",  help="Use all data available provided: kma databases, external db if provided, genbank, user data..")
exclusive_group2.add_argument("--only_user_data", action="store_true",  help="Use only indexed genomes previously identified stored in database folder.")
exclusive_group2.add_argument("--only_genbank_data", action="store_true",  help="Use only reference genomes previously downloaded from NCBI stored in database folder.")
exclusive_group2.add_argument("--only_external_kma", action="store_true", help="Provide this option if only database(s) provided via --kma_external_file should be used for species identification.")

parameters_group_ident = subparser_ident.add_argument_group("Parameters")
parameters_group_ident.add_argument("--threads", type=int, help="Number of CPUs to use [Default: 2].", default=2)
parameters_group_ident.add_argument("--KMA_cutoff", type=int, help="Similarity cutoff for databases filtering. Range: 1-100. [Default = 80]", default=80)
parameters_group_ident.add_argument("--fast", action="store_true", help="Do not update database, just identify samples according to databases provided [Default OFF].")

info_group_ident = subparser_ident.add_argument_group("Additional information")
info_group_ident.add_argument("--debug", action="store_true", help="Show additional message for debugging purposes.")
info_group_ident.add_argument("--help_format", action="store_true", help="Show additional help on name format for files.")
info_group_ident.add_argument("--help_project", action="store_true", help="Show additional help on the project scheme.")
info_group_ident.add_argument("--help_KMA", action="store_true", help="Show additional help on KMA software and options.")
info_group_ident.add_argument("--help_MLSTar", action="store_true", help="Show additional help on MLSTar software and options.")

subparser_ident.set_defaults(func=BacterialTyper.modules.ident.run_ident)
##-------------------------------------------------------------##

##--------------------------- profile ------------------------ ##
subparser_profile = subparsers.add_parser(
    'profile',
    help='Virulence & Resistance profile.',
    description='This module generates a resistance and virulence profile using several databases ...',
)
profile_mode_name = subparser_profile.add_argument_group("Mode")
profile_mode = profile_mode_name.add_mutually_exclusive_group(required= not any(elem in help_options for elem in sys.argv))
profile_mode.add_argument("--project", action="store_true", help="Requires as --input a folder containing a project with samples, metadata, configuration etc. [Default]")
profile_mode.add_argument("--detached", action="store_true", help="Isolated mode. --input is a folder containing assemblies or scaffolds. Provide a unique path o several using --batch option")

initial_group_profile = subparser_profile.add_argument_group("Input/Output")
initial_group_profile.add_argument("--input", help="Folder containing trimmed reads. Files could be .fastq/.fq/ or fastq.gz/.fq.gz", required= not any(elem in help_options for elem in sys.argv))
initial_group_profile.add_argument("--output_folder", help="Output folder.", required = '--detached' in sys.argv)
initial_group_profile.add_argument("--single_end", action="store_true", help="Single end files [Default OFF]. Default mode is paired-end.")
initial_group_profile.add_argument("--batch", action="store_true", help="Provide this option if input is a file containing multiple paths instead a path.")
initial_group_profile.add_argument("--in_sample", help="File containing a list of samples to include (one per line) from input folder(s) [Default OFF].")
initial_group_profile.add_argument("--ex_sample", help="File containing a list of samples to exclude (one per line) from input folder(s) [Default OFF].")

dataset_group_profile = subparser_profile.add_argument_group("Datasets")
dataset_group_profile.add_argument("--database", help="Directory containing databases previously downloaded such as ARIBA, KMA, BUSCO genbank and user_data folders.", required=not any(elem in help_options for elem in sys.argv))
dataset_group_profile.add_argument("--additional_database", action="store_true", help="## [TODO]")

parameters_group_profile = subparser_profile.add_argument_group("Parameters")
parameters_group_profile.add_argument("--threads", type=int, help="Number of CPUs to use [Default: 2].", default=2)
parameters_group_profile.add_argument("--fast", action="store_true", help="Do not update database, just identify samples according to databases provided [Default OFF].")

## ariba
ariba_db_group_profile = subparser_profile.add_argument_group("ARIBA databases")
ariba_db_group_profile.add_argument("--ARIBA_db", dest='ariba_dbs', nargs='*', help="ARIBA database(s) to download. Provide several input if desired. Not compatible with --no_ARIBA. [Default: CARD & VFDB]", choices=['CARD', 'MEGARes', 'VFDB', 'ResFinder', 'PlasmidFinder', 'VirulenceFinder', 'srst2', 'ARG-ANNOT'])
ariba_db_group_profile.add_argument("--no_def_ARIBA", action="store_true", help="Only applicable if ARIBA_db is ON. It discards default databases (CARD & VFDB) from ARIBA_db. Only user selected databases indexed. [Default: OFF]")
ariba_db_group_profile.add_argument("--ARIBA_cutoff", type=float, help="ARIBA assembly threshold cutoff [0-1]. [ Default: 0.90]", default=0.90)

info_group_profile = subparser_profile.add_argument_group("Additional information")
info_group_profile.add_argument("--debug", action="store_true", help="Show additional message for debugging purposes.")
info_group_profile.add_argument("--help_format", action="store_true", help="Show additional help on name format for files.")
info_group_profile.add_argument("--help_project", action="store_true", help="Show additional help on the project scheme.")
info_group_profile.add_argument("--help_ARIBA", action="store_true", help="Print further information for ARIBA databases")

subparser_profile.set_defaults(func=BacterialTyper.modules.profile.run_profile)

##-------------------------------------------------------------##

##--------------------------- MGE ---------------------------- ##
subparser_MGE = subparsers.add_parser(
    'MGE',
    help='Mobile Genetic Elements (MGE) analysis.',
    description='This module identifies Mobile Genetic elements: plasmids, bacteriophages and genomic islands',
)
MGE_mode_name = subparser_MGE.add_argument_group("Mode")
MGE_mode = MGE_mode_name.add_mutually_exclusive_group(required= not any(elem in help_options for elem in sys.argv))
MGE_mode.add_argument("--project", action="store_true", help="Requires as --input a folder containing a project with samples, metadata, configuration etc. [Default]")
MGE_mode.add_argument("--detached", action="store_true", help="Isolated mode. Provided as --input a csv file containing different fields. See --help_input_MGE for further details.")

initial_group_MGE = subparser_MGE.add_argument_group("Input/Output")
initial_group_MGE.add_argument("--input", help="Folder containing assemblies. Contig/Scaffolds files could be within a project folder or in the input folder. File must end with tag '_chromosome(.fasta/.fna)'.\nIf not, provide full path using --batch option", required = '--project' in sys.argv)
initial_group_MGE.add_argument("--output_folder", help="Output folder.", required = '--detached' in sys.argv)
initial_group_MGE.add_argument("--single_end", action="store_true", help="Single end files [Default OFF]. Default mode is paired-end.")
initial_group_MGE.add_argument("--in_sample", help="File containing a list of samples to include (one per line) from input folder(s) [Default OFF].")
initial_group_MGE.add_argument("--ex_sample", help="File containing a list of samples to exclude (one per line) from input folder(s) [Default OFF].")

parameters_group_MGE = subparser_MGE.add_argument_group("Parameters")
parameters_group_MGE.add_argument("--threads", type=int, help="Number of CPUs to use [Default: 2].", default=2)
parameters_group_MGE.add_argument("--database", help="Directory containing databases previously downloaded such as ARIBA, KMA, BUSCO genbank and user_data folders.", required= not any(elem in help_options for elem in sys.argv))

options_group_MGE_name = subparser_MGE.add_argument_group("Analysis")
options_group_MGE = options_group_MGE_name.add_mutually_exclusive_group(required= not any(elem in help_options for elem in sys.argv))
options_group_MGE.add_argument("--plasmid_analysis", action="store_true", help="Identify putative plasmids")
options_group_MGE.add_argument("--phage_analysis", action="store_true", help="Identify putative phages inserted in the genome")
options_group_MGE.add_argument("--GI_analysis", action="store_true", help="Identify putative Genomic Islands.")
options_group_MGE.add_argument("--all_data", action="store_true", help="Identify GIs, phages and plasmids.")

phispy_options = subparser_MGE.add_argument_group("PhiSpy options")
phispy_options.add_argument("--training_set", type=int, help="Choose a training set from the available training sets [Default: 0 (Generic)]. See additional details with --help_PhiSpy", default=0)
phispy_options.add_argument("--window_size", type=int, help="Window size of consecutive genes to look through to find phages. [Default: 20]", default=20)
phispy_options.add_argument("--phage_genes", type=int, help="Number of consecutive genes in a region of window size that must be prophage genes to be called. [Default: 5]", default=5)

info_group_MGE = subparser_MGE.add_argument_group("Additional information")
info_group_MGE.add_argument("--debug", action="store_true", help="Show additional message for debugging purposes.")
info_group_MGE.add_argument("--help_format", action="store_true", help="Show additional help on name format for files.")
info_group_MGE.add_argument("--help_project", action="store_true", help="Show additional help on the project scheme.")
info_group_MGE.add_argument("--help_input_MGE", action="store_true", help="Print further information for input options under MGE module.")
info_group_MGE.add_argument("--help_MGE_analysis", action="store_true", help="Print further information for Mobile Genetic Element module analysis.")
info_group_MGE.add_argument("--help_PhiSpy", action="store_true", help="Print further information for PhiSpy analysis.")

subparser_MGE.set_defaults(func=BacterialTyper.modules.MGE.run_MGE)
##-------------------------------------------------------------##

##--------------------------- cluster ---------------------------- ##
subparser_cluster = subparsers.add_parser(
    'cluster',
    help='Cluster sequence analysis.',
    description='This module calls ...',
)
cluster_mode_name = subparser_cluster.add_argument_group("Mode")
cluster_mode = cluster_mode_name.add_mutually_exclusive_group(required= not any(elem in help_options for elem in sys.argv))
cluster_mode.add_argument("--project", action="store_true", help="Requires as --input a folder containing a project with samples, metadata, configuration etc. [Default]")
cluster_mode.add_argument("--detached", action="store_true", help="Isolated mode. Provided as --input a csv file containing different fields. See --help_input_MGE for further details.")

initial_group_cluster = subparser_cluster.add_argument_group("Input/Output")
initial_group_cluster.add_argument("--input", help="Folder containing assemblies. Contig/Scaffolds files could be within a project folder or in the input folder. File must end with tag '_chromosome(.fasta/.fna)'.\nIf not, provide full path using --batch option", required = '--project' in sys.argv)
initial_group_cluster.add_argument("--output_folder", help="Output folder.", required = '--detached' in sys.argv)
initial_group_cluster.add_argument("--in_sample", help="File containing a list of samples to include (one per line) from input folder(s) [Default OFF].")
initial_group_cluster.add_argument("--ex_sample", help="File containing a list of samples to exclude (one per line) from input folder(s) [Default OFF].")
initial_group_cluster.add_argument("--batch", action="store_true", help="Provide this option if input is a file containing multiple paths instead a path.")

parameters_group_cluster = subparser_cluster.add_argument_group("Parameters")
parameters_group_cluster.add_argument("--threads", type=int, help="Number of CPUs to use [Default: 2].", default=2)
parameters_group_cluster.add_argument("--database", help="Directory containing databases previously downloaded such as ARIBA, KMA, BUSCO genbank and user_data folders.", required= not any(elem in help_options for elem in sys.argv))
parameters_group_cluster.add_argument("--external_file", help="Fasta file to include in the clustering.")
parameters_group_cluster.add_argument("--batch_external", action="store_true", help="Provide this option if file provided via --external_file is a file containing multiple paths instead a path.")

options_group_cluster_name = subparser_cluster.add_argument_group("Analysis")
options_group_cluster = options_group_cluster_name.add_mutually_exclusive_group(required= not any(elem in help_options for elem in sys.argv))
options_group_cluster.add_argument("--all_data", action="store_true",  help="Use all data available provided: project samples, genbank, previous identified user data..")
options_group_cluster.add_argument("--only_project_data", action="store_true", help="Only use samples in the project folder provided.")
options_group_cluster.add_argument("--user_data", action="store_true",  help="Use only indexed genomes previously identified stored in database folder.")
options_group_cluster.add_argument("--genbank_data", action="store_true",  help="Use only reference genomes previously downloaded from NCBI stored in database folder.")
options_group_cluster.add_argument("--only_external_data", action="store_true",  help="Use only reference genomes previously downloaded from NCBI stored in database folder.")

parameters_minHash_group_cluster = subparser_cluster.add_argument_group("Parameters MinHash")
parameters_minHash_group_cluster.add_argument("--n_sketch", type=int, help="Sketch size. Each sketch will have at most this many non-redundant min-hashes. [Default: 5000].", default=5000)
parameters_minHash_group_cluster.add_argument("--kmer_size", type=int, help="Hashes will be based on strings of this many nucleotides. [Default: 51]", default=51)

info_group_cluster = subparser_cluster.add_argument_group("Additional information")
info_group_cluster.add_argument("--debug", action="store_true", help="Show additional message for debugging purposes.")
info_group_cluster.add_argument("--help_format", action="store_true", help="Show additional help on name format for files.")
info_group_cluster.add_argument("--help_project", action="store_true", help="Show additional help on the project scheme.")
info_group_cluster.add_argument("--help_Mash", action="store_true", help="Print further information for Mash analysis.")

subparser_cluster.set_defaults(func=BacterialTyper.modules.cluster.run_cluster)
##-------------------------------------------------------------##

##--------------------------- cluster ---------------------------- ##
subparser_phylo = subparsers.add_parser(
    'phylo',
    help='Phylogenetic analysis.',
    description='This module calls ...',
)
phylo_mode_name = subparser_phylo.add_argument_group("Mode")
phylo_mode = phylo_mode_name.add_mutually_exclusive_group(required= not any(elem in help_options for elem in sys.argv))
phylo_mode.add_argument("--project", action="store_true", help="Requires as --input a folder containing a project with samples, metadata, configuration etc. [Default]")
phylo_mode.add_argument("--detached", action="store_true", help="Isolated mode. Provided as --input a csv file containing different fields. See --help_input_MGE for further details.")

initial_group_phylo = subparser_phylo.add_argument_group("Input/Output")
initial_group_phylo.add_argument("--input", help="Folder containing assemblies. Contig/Scaffolds files could be within a project folder or in the input folder. File must end with tag '_chromosome(.fasta/.fna)'.\nIf not, provide full path using --batch option", required = '--project' in sys.argv)
initial_group_phylo.add_argument("--output_folder", help="Output folder.", required = '--detached' in sys.argv)
initial_group_phylo.add_argument("--in_sample", help="File containing a list of samples to include (one per line) from input folder(s) [Default OFF].")
initial_group_phylo.add_argument("--ex_sample", help="File containing a list of samples to exclude (one per line) from input folder(s) [Default OFF].")
initial_group_phylo.add_argument("--batch", action="store_true", help="Provide this option if input is a file containing multiple paths instead a path.")
initial_group_phylo.add_argument("--single_end", action="store_true", help="Single end files [Default OFF]. Default mode is paired-end.")

reference_group_phylo_name = subparser_phylo.add_argument_group("Reference")
reference_group_phylo = reference_group_phylo_name.add_mutually_exclusive_group(required= not any(elem in help_options for elem in sys.argv))
reference_group_phylo.add_argument("--Genbank_ID", help="Genbank ID of the reference strain. Available or not in database")
reference_group_phylo.add_argument("--user_sample_ID", help="Sample ID previously analyzed and available in database.")
reference_group_phylo.add_argument("--project_sample_ID", help="Project sample ID.")
reference_group_phylo.add_argument("--user_gbk", help="Genbank file format provided by the user.")

parameters_group_phylo = subparser_phylo.add_argument_group("Parameters")
parameters_group_phylo.add_argument("--threads", type=int, help="Number of CPUs to use [Default: 2].", default=2)
parameters_group_phylo.add_argument("--database", help="Directory containing databases previously downloaded such as ARIBA, KMA, BUSCO genbank and user_data folders.", required= not any(elem in help_options for elem in sys.argv))
parameters_group_phylo.add_argument("--name", help="Name ID to identify the analysis", required= not any(elem in help_options for elem in sys.argv))
parameters_group_phylo.add_argument("--other_options", help="String of options to include in snippy call")
parameters_group_phylo.add_argument("--output_format", help="Alignment output format. [Default: fasta]", default="fasta", choices=['nexus', 'phylip', 'clustalw', 'fasta'])

options_group_phylo_name = subparser_phylo.add_argument_group("Analysis")
options_group_phylo = options_group_phylo_name.add_mutually_exclusive_group(required= not any(elem in help_options for elem in sys.argv))
options_group_phylo.add_argument("--all_data", action="store_true",  help="Use all data available provided: project samples, genbank, previous identified user data..")
options_group_phylo.add_argument("--only_project_data", action="store_true", help="Only use samples in the project folder provided.")
options_group_phylo.add_argument("--user_data", action="store_true",  help="Use only indexed genomes previously identified stored in database folder.")
options_group_phylo.add_argument("--genbank_data", action="store_true",  help="Use only reference genomes previously downloaded from NCBI stored in database folder.")

info_group_phylo = subparser_phylo.add_argument_group("Additional information")
info_group_phylo.add_argument("--debug", action="store_true", help="Show additional message for debugging purposes.")
info_group_phylo.add_argument("--help_format", action="store_true", help="Show additional help on name format for files.")
info_group_phylo.add_argument("--help_project", action="store_true", help="Show additional help on the project scheme.")
info_group_phylo.add_argument("--help_Snippy", action="store_true", help="Print further information for Snippy analysis.")

subparser_phylo.set_defaults(func=BacterialTyper.modules.phylo.run_phylo)
##-------------------------------------------------------------##




## space
subparser_space = subparsers.add_parser(' ', help='')

#############################
#### Additional messages ####
#############################
##------------------------------ Information help ---------------------- ##
subparser_help = subparsers.add_parser(
    'info',
    help='Print additional information & help messages ',
    description='For different modules, options or parameters print additional information and help messages',
)
subparser_help_name = subparser_help.add_argument_group("Show additional help information")
subparser_help = subparser_help_name.add_mutually_exclusive_group(required= True)
subparser_help.add_argument("--help_project", action="store_true", help="...")
subparser_help.add_argument("--help_format", action="store_true", help="...")
subparser_help.add_argument("--help_trimm_adapters", action="store_true", help="Show additional information on trimm adapters.")
subparser_help.add_argument("--help_BUSCO", action="store_true", help="...")
subparser_help.add_argument("--help_Prokka", action="store_true", help="...")
subparser_help.add_argument("--help_multiqc", action="store_true", help="Show additional help on the multiQC module.")
subparser_help.add_argument("--help_ARIBA", action="store_true", help="...")
subparser_help.add_argument("--help_KMA", action="store_true", help="Show additional help on KMA software and options.")
subparser_help.add_argument("--help_input_MGE", action="store_true", help="Print further information for input options under MGE module.")
subparser_help.add_argument("--help_MGE_analysis", action="store_true", help="Print further information for Mobile Genetic Element module analysis.")
subparser_help.add_argument("--help_PhiSpy", action="store_true", help="Print further information for PhiSpy analysis.")
subparser_help.add_argument("--help_Mash", action="store_true", help="Print further information for Mash clustering analysis.")
subparser_help.add_argument("--help_Snippy", action="store_true", help="Print further information for Snippy analysis.")
subparser_help.add_argument("--help_MLSTar", action="store_true", help="Print further information for MLST analysis.")

subparser_help.set_defaults(func=BacterialTyper.modules.help_info.run_info)
##-------------------------------------------------------------##

##-------------------------- version ------------------------- ##
subparser_version = subparsers.add_parser(
    'version',
    help='Packages & software versions.',
    description='This code generates prints an index of version number for the different packages and other softwares employed here',
)
subparser_version.add_argument("option", help="Print only this pipeline version or all packages versions.", choices=['only','all'])
subparser_version.set_defaults(func=BacterialTyper.modules.version.run)
##-------------------------------------------------------------##

##--------------------------- citation ------------------------##
subparser_citation = subparsers.add_parser(
    'citation',
    help='Packages & software citations.',
    description='This code prints an index of citation for the different packages and other softwares employed here',
)
subparser_citation.add_argument("option", help="Print only this pipeline citation or all packages references.", choices=['only','all'])
subparser_citation.set_defaults(func=BacterialTyper.modules.citation.run)
##-------------------------------------------------------------##

subparser_space = subparsers.add_parser(' ', help='')

#####
args = parser.parse_args()
if hasattr(args, 'func'):
    args.func(args)
else:
    parser.print_help()
