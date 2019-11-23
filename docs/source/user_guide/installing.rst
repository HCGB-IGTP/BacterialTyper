.. _installing:

Installation
************

.. contents::

This is an installation guide for BacterialTyper. 

First get the code from different sources available and then 
make sure you have all the dependencies fulfilled from section :ref:`Requirements-dependencies`.

We encourage you to create a virtual environment before installing BacterialTyper.

.. note::

    If you wish to contribute to the project, it's recommended you
       :ref:`install the latest development version<install-from-source>`.

.. _virtual-env-BacterialTyper:

Python environment
==================

It is recommended to install BacterialTyper and all its dependencies and modules under a python environment. 

.. ################################################
.. create-env:

Create environment
------------------
To do so, run :file:`config/main/python_environment-commands.sh`, which is shown below:

.. include:: ../../../config/main/python_environment-commands.sh
   :literal:
   
Execute the file from the directory :file:`BacterialTyper/` as::

   sh config/main/python_environment-commands.sh
   
.. ################################################
.. _activate-env:

Activate environemt
-------------------

Before executing any BacterialTyper module or script, it is convenient to activate the environment.

To do so, run :file:`config/main/activate_environment.sh`, which is shown below:

.. include:: ../../../config/main/activate_environment.sh
   :literal:

Execute the file from the directory :file:`BacterialTyper/` as::
   sh config/main/activate_environment.sh

.. ################################################
.. _install-BacterialTyper:

Installing BacterialTyper
=========================

If you want to run a stable version of BacterialTyper, install it using pip package and following instructions
in :ref:`install-from-pip`.

If you are interested in contributing to BacterialTyper development, running the latest source code, 
or just like to build everything yourself, it is not difficult to build it from source following 
the :ref:`install-from-source` instructions. 

.. ################################################
.. _Requirements-dependencies:

Requirements and dependencies
-----------------------------

There are several extra python module requirements that are needed. 
They are listed in :file:`config/main/python_requirements.txt`, which is shown below:

.. include:: ../../../config/main/python_requirements.txt
   :literal:
   
Additionally, several third-party software packages are also required. 
They are listed in :file:`config/main/DEPENDENCIES`, which is shown below:

.. include:: ../../../config/main/DEPENDENCIES
   :literal:

.. ################################################
.. _install-from-pip:

Installing from pip
-------------------

Follow the `pip installing instructions`_ to learn about installing packages.

Type::
   python -m pip install BacterialTyper

.. ################################################
.. _install-from-source:

Installing from source
----------------------

To get the code, you can either grab the latest tar.gz release file from the `PyPI files page`_ and follow the instructions 
in :ref:`get-Pypi-targz` and install it. On the other hand, if you want to develop BacterialTyper or just 
need the latest bugfixed version, grab the latest git version from the `BacterialTyper github`_ website, and 
follow the :ref:`get-git-code` instructions.

.. _get-Pypi-targz: 

Get latest Pypi files
^^^^^^^^^^^^^^^^^^^^^
   
The Python Package Index (PyPI) is a repository of software for the Python programming language. PyPI helps 
you find and install software developed and shared by the Python community. 

TODO: Fill it.

.. _get-git-code: 

Get latest git version
^^^^^^^^^^^^^^^^^^^^^^

Type:
   git clone https://github.com/JFsanchezherrero/BacterialTyper.git 

Add git clone submodules command

.. include:: ../links.inc
