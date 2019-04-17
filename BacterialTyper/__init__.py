#from pkg_resources import get_distribution
#try:
#    __version__ = get_distribution('BacterialTyper').version
#except:
#    __version__ = 'local'
### to include when distribution available
## version will be retrieve from setup.py


__all__ = [
	'config',
	'modules',
	'annotation',
	'bacteriophage',
	'functions',
	'blast_parser',
	'sampleParser',
	'config',
	'MLSTar',
	'trimmomatic_call',
	'database_generator',
	'fastqc_caller',
	'spades_assembler',
	'species_identification_KMA',
	'extern_progs',
	'multiQC_report'
	
]

from BacterialTyper import *




