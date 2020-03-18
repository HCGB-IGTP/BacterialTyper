import os
import shutil
import sys
import glob
from setuptools import setup, find_packages, Extension

def get_require_modules():
    """
    Get main python requirements modules
    """
    with open("./docs/source/data/python_requirements_summary.csv", 'r') as f:
        myModules = [line.strip() for line in f]
    
    ## To Do: Debug
    print(myModules)
    return myModules

def get_version():
    """
    Original code: PhiSpy setup.py 
    https://github.com/linsalrob/PhiSpy/blob/master/setup.py
    """
    with open("VERSION", 'r') as f:
        v = f.readline().strip()
    return v

### subdmodules:: Islandpath-DIMOB
islandpath_files = [os.path.join('third_party', 'Islandpath-DIMOB', x) for x in islandpath_files]
islandpath_mod = Extension(
    "Islandpath-DIMOB",
    islandpath_files,
    include_dirs=[os.path.join('third_party', 'Islandpath-DIMOB')],
)

setup(
	name='BacterialTyper',
    version=get_version(),
    description='xxxxx',
    packages = find_packages(),
    author='Jose F Sanchez Herrero',
    author_email='xxxx',
    url='xxxx',
    scripts=["main/BacterialTyper.py"],
    install_requires=get_require_modules(),
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
