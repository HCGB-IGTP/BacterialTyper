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
    'card_trick_caller',
    'database_generator',
    'database_user',
    'edirect_caller',
    'fastqc_caller',
    'functions',
    'genomic_island',
    'MLSTar',
    'min_hash_caller',
    'multiQC_report',
    'sampleParser',
    'spades_assembler',
    'species_identification_KMA',
    'trimmomatic_call',
    'variant_calling',
    'virulence_resistance'
]

from BacterialTyper.scripts import *