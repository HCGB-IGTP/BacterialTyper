#from pkg_resources import get_distribution
#try:
#    __version__ = get_distribution('BacterialTyper').version
#except:
#    __version__ = 'local'
### to include when distribution available
## version will be retrieve from setup.py


__all__ = [
	'BUSCO_caller',
	'MLSTar',
	'annotation',
	'ariba_caller',
	'bacteriophage',
	'blast_parser',
	'config',
	'database_generator',
	'extern_progs',
	'fastqc_caller',
	'functions',
	'modules',
	'multiQC_report',
	'other_tools',
	'sampleParser',
	'set_configuration',
	'spades_assembler',
	'species_identification_KMA',
	'trimmomatic_call',
	'virulence_resistance'
]

from BacterialTyper import *




