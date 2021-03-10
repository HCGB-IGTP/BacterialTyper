# Create pip package


## clean distribution packages
'''
sh devel/pypi/clean_devel.sh
'''

## create distribution files
'''
sh devel/pypi/create_distro.sh
'''

## create .pypirc file

'''
$ nano .pypirc
[distutils] 
index-servers=pypi
[pypi] 
repository = https://upload.pypi.org/legacy/ 
username =jfsanchezherrero
'''

## Upload using twine
'''
sh devel/pypi/upload_pypi.sh
'''

## references
https://dzone.com/articles/executable-package-pip-install
https://packaging.python.org/tutorials/packaging-projects/
conda build -c bioconda conda
