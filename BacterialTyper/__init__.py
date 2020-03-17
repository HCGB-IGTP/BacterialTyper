#from pkg_resources import get_distribution
#try:
#    __version__ = get_distribution('BacterialTyper').version
#except:
#    __version__ = 'local'
### to include when distribution available
## version will be retrieve from setup.py

__all__ = [
	'data',
	'modules',
	'other_tools',
	'scripts'
]

from BacterialTyper import *




