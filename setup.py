import os
import shutil
import sys
import glob
import setuptools

#######
def get_require_modules():
    """
    Get main python requirements modules
    """
    with open("./BacterialTyper/config/python/python_requirements_summary.csv", 'r') as f:
        myModules = [line.strip().split(',')[0] for line in f]
    
    return myModules

#######
def get_version():
    """
    Original code: PhiSpy setup.py 
    https://github.com/linsalrob/PhiSpy/blob/master/setup.py
    """
    with open("./VERSION", 'r') as f:
        v = f.readline().strip()
    return v

### subdmodules:: Islandpath-DIMOB
#islandpath_files = [os.path.join('third_party', 'Islandpath-DIMOB', x) for x in islandpath_files]
#islandpath_mod = Extension(
#    "Islandpath-DIMOB",
#    islandpath_files,
#    include_dirs=[os.path.join('third_party', 'Islandpath-DIMOB')],
#)
#######
long_description_text = ""
with open("README.md", "r") as fh:
    long_description_text = fh.read()

#######
setuptools.setup(
    name='BacterialTyper',
    version=get_version(),

    scripts=glob.glob('main/*'),
    author="Jose F. Sanchez-Herrero",
    author_email="jfbioinformatics@gmail.com",

    long_description_content_type="text/markdown",
    long_description=long_description_text,

    url="https://github.com/HCGB-IGTP/BacterialTyper/",

    install_requires=get_require_modules(),
    packages = setuptools.find_packages(),
    include_package_data=True,

    license='MIT License',    
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent"
    ],
)

## https://github.com/sanger-pathogens/ariba/blob/master/setup.py
## https://dzone.com/articles/executable-package-pip-install
## https://github.com/BillMills/pythonPackageLesson
