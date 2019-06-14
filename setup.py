import os
import shutil
import sys
import glob
from setuptools import setup, find_packages, Extension


PhiSpy_files = [os.path.join('third_party', 'PhiSpy', x) for x in PhiSpy_files]
PhiSpy_mod = Extension(
    "PhiSpy",
    PhiSpy_files,
    include_dirs=[os.path.join('third_party', 'PhiSpy')],
)

setup(
	ext_modules=[PhiSpy_mod],
    name='xxxx',
    version='xx.xx.xx',
    description='xxxxx',
    packages = find_packages(),
    author='Jose F Sanchez Herrero',
    author_email='xxxx',
    url='xxxx',
    scripts=glob.glob('scripts/*'),
    install_requires=[
		## add additional modules
        'BeautifulSoup4 >= 4.1.0', 
        'biopython',
        'dendropy >= 4.2.0',
        'matplotlib',
        'pyfastaq >= 3.12.0',
        'pysam >= 0.9.1',
        'pymummer<=0.10.3',
    ],
    license='GPLv3',    
    classifiers=[
        'Development Status :: 4 - Beta',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Programming Language :: Python :: 3 :: Only',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
    ],
)

## https://github.com/sanger-pathogens/ariba/blob/master/setup.py
## https://dzone.com/articles/executable-package-pip-install
## https://github.com/BillMills/pythonPackageLesson
