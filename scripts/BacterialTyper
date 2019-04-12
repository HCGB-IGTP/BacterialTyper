 #usr/bin/env python
'''
This code...
Jose F. Sanchez
Copyright (C) 2019 Lauro Sumoy Lab, IGTP, Spain
'''
import argparse
import os
import BacterialTyper

## useful imports
#pythonDir = os.path.dirname(os.path.realpath(__file__)) + '/tools/python/'

parser = argparse.ArgumentParser(
    prog='BacterialGenotyper',
    usage='BacterialGenotyper <command> <options>',
    description='BacterialGenotyper: Bacterial Genotyping using NGS...',

)
subparsers = parser.add_subparsers(title='Available commands', help='', metavar='')

subparser_database = subparsers.add_parser(
    'database_generator',
    help='Downloads, updates and prepares database for later use.',
    usage='BacterialGenotyper database_generator <option> <folder> [ID_file]',
    description='....',
)
#database_py = pythonDir + '/database_generator.py'
subparser_database.set_defaults(func=BacterialTyper.scripts.database_generator.help_options)

subparser_aln2meta = subparsers.add_parser(
    'aln2meta',
    help='Converts multi-aln fasta and SNPs to metadata',
    usage='ariba aln2meta [options] <aln_fasta> <variants_tsv> <(non)coding> <outprefix>',
    description='Make metadata input to prepareref, using multialignment and SNPs',
)

args = parser.parse_args()

if hasattr(args, 'func'):
    args.func(args)
else:
    parser.print_help()
