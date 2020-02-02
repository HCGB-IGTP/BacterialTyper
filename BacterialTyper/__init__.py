#from pkg_resources import get_distribution
#try:
#    __version__ = get_distribution('BacterialTyper').version
#except:
#    __version__ = 'local'
### to include when distribution available
## version will be retrieve from setup.py


__all__ = [
	'annotation',
	'ariba_caller',
	'bacteriophage',
	'blast_parser',
	'BUSCO_caller',
	'config',
	'card_trick_caller',
	'data',
	'database_generator',
	'database_user',
	'edirect_caller',
	'extern_progs',
	'fastqc_caller',
	'functions',
	'genomic_island',
	'MLSTar',
	'min_hash_caller',
	'modules',
	'multiQC_report',
	'other_tools',
	'sampleParser',
	'spades_assembler',
	'species_identification_KMA',
	'trimmomatic_call',
	'variant_calling',
	'virulence_resistance'
]

from BacterialTyper import *




