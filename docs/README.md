# README Documentation

This README contains information regarding the documentation for BacterialTyper

The BacterialTyper documentation is supported and deployed using Sphinx and rinohtype.

Instruction and examples on how to create this documentation were obtained from:
- https://github.com/richdayandnight/Tutorial_SimpleTeacherAPI
- https://pythonhosted.org/an_example_pypi_project/sphinx.html#is-sweaty
- https://medium.com/@richdayandnight/a-simple-tutorial-on-how-to-document-your-python-project-using-sphinx-and-rinohtype-177c22a15b5b
- http://www.sphinx-doc.org/en/master/
- http://www.sphinx-doc.org/en/master/usage/quickstart.html
- https://github.com/matplotlib/sampledoc

See configuration required in `config/` folder for python, perl and system packages to be installed.

- Python
```
sphinx
sphinxmark
rinohtype
sphinxcontrib-bibtex
```

- Perl
```
latexmk,http://mirrors.ctan.org/support/latexmk.zip
```

- System packages
```
sudo apt install texlive-full
sudo apt install texmaker
sudo apt install graphviz
```


To create the documentation execute:

```
sphinx-quickstart
make html
sphinx-build -b rinoh source build/rinoh
```
