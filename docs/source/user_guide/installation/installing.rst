.. ################################################
.. _installing:

Installation
************

.. toctree::
   :hidden:
   :maxdepth: 1
   
   requirements.rst

.. contents::

This is an installation guide for BacterialTyper. 

First get the code from different sources available and then 
make sure you have all the dependencies fulfilled from section :ref:`Requirements-dependencies`.

.. note::

    If you want to contribute or see additional details of the project, it's recommended you
    :ref:`install the latest development version<install-from-source>`.


.. ################################################
.. _Requirements-dependencies:
.. ################################################

Requirements and dependencies
=============================

BacterialTyper is a pipeline composed of multiple lines of code that call and automatizes the analysis from
multiple bioinformatic tools. Several dependencies (python, perl, third party software) are required for
the BacterialTyper analysis. 

We include details on the different modules required here:

.. toctree:: 
   requirements.rst

.. ################################################
.. _install-BacterialTyper:
.. ################################################

Installing BacterialTyper
=========================

If you want to run a stable version of BacterialTyper, install it using pip package and following instructions
in :ref:`install-from-pip`.

If you are interested in contributing to BacterialTyper development, running the latest source code, 
or just like to build everything yourself, it is not difficult to build it from source following 
the :ref:`install-from-source` instructions. 

  
.. ##################
.. _install-from-pip:
.. ##################

Installing from pip
-------------------

The fastest and most straightforward solution to install the latest release (or any desired version) of BacterialTyper 
is to install it and all its dependendencies, using `Python pip`_. Follow the `pip installing instructions`_ to 
learn about installing packages.

Type:

.. code-block:: sh

   python -m pip install BacterialTyper


First, we encourage you to create a virtual environment before installing BacterialTyper. See 
section :ref:`virtual-env-BacterialTyper` for further details.


.. ##################
.. _install-from-source:
.. ##################

Installing from source
----------------------

To get the code, you can either grab the latest tar.gz release file from the `PyPI files page`_ and follow the instructions 
in :ref:`get-Pypi-targz` and install it. On the other hand, if you want to develop BacterialTyper or just 
need the latest bugfixed version, grab the latest git version from the `BacterialTyper github`_ website, and 
follow the :ref:`get-git-code` instructions.


Get source code
^^^^^^^^^^^^^^^


.. _get-Pypi-targz: 

Get latest Pypi files
"""""""""""""""""""""
   
The Python Package Index (PyPI) is a repository of software for the Python programming language. PyPI helps 
you find and install software developed and shared by the Python community. 

TODO: Fill it.

.. _get-git-code: 

Get latest git version
""""""""""""""""""""""

Type:
   git clone https://github.com/JFsanchezherrero/BacterialTyper.git 

Add git clone submodules command

Install Dependencies
^^^^^^^^^^^^^^^^^^^^

.. ################################################
.. _install-perl_packages:

Install perl packages
"""""""""""""""""""""

There is shell script available for the perl package modules installation (file :file:`config/main/perl_lib_installer.sh`). 

.. include:: ../../../../config/main/perl_lib_installer.sh
   :literal:
 
Execute the file from the directory :file:`BacterialTyper/` as:

.. code-block:: sh

   sh config/main/perl_lib_installer.sh



.. ############################
.. _virtual-env-BacterialTyper:
.. ############################

Python environment
==================

It is recommended to install BacterialTyper and all its dependencies and modules under a python environment. 

.. ###########
.. create-env:
.. ###########

Create environment
------------------
To do so, run :file:`config/main/python_environment-commands.sh`, which is shown below:

.. include:: ../../../../config/main/python_environment-commands.sh
   :literal:
   
Execute the file from the directory :file:`BacterialTyper/` as:

.. code-block:: sh

   sh config/main/python_environment-commands.sh
   
.. ###########
.. _activate-env:
.. ###########

Activate environemt
-------------------

Before executing any BacterialTyper module or script, it is convenient to activate the environment.

To do so, run :file:`config/main/activate_environment.sh`, which is shown below:

.. include:: ../../../../config/main/activate_environment.sh
   :literal:

Execute the file from the directory :file:`BacterialTyper/` as:

.. code-block:: sh
  
   sh config/main/activate_environment.sh




.. include:: ../../links.inc

