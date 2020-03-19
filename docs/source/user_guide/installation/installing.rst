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

Under some circumstancies (develop, bugs fixed, etc) you might be interested in obtaining the latest code
 version. Take into account, that you will to install dependencies and fulfill requirements to have a working
 distribution. 

.. ##############
.. _get-git-code:
.. ##############

Get source code
^^^^^^^^^^^^^^^

The BacterialTyper project uses git_ as a version control system. To get the code, you can grab the latest version 
from the `BacterialTyper github`_ website, and follow the :ref:`install-BacterialTyper-source` instructions.

Using the command-line, check you have a working distribution of git by typing ``git --help`` or install it by typing:

.. code-block:: sh

   sudo apt update
   sudo apt upgrade
   sudo apt install git

Once you have ``git`` installed and working change directory to a folder of interest where the BacterialTyper 
project would be download. Type:

.. code-block:: sh

   git clone https://github.com/JFsanchezherrero/BacterialTyper.git 

Additionally, BacterialTyper contains a third-party project, IslandPath-DIMOB (IslandPath_ :cite:`Bertelli2017`) that BacterialTyper developers
forked, updated and adapted. See additional details in section: :ref:`islandPath-DIMOB-forked`.

Then, if you are downloading BacterialTyper code from source, you need to retrieve the submodule latest version
and install it accordingly. To do so, type:

.. code-block:: sh
   
   git clone --recurse-submodules https://github.com/JFsanchezherrero/BacterialTyper.git

On the other hand, you can also do it manually. In the directory containing git BacterialTyper project you would 
find a ``third_party`` folder. Change directory to there and type:

.. code-block:: sh

   cd BacterialTyper/third_party
   git submodule init
   git submodule update
   
For additional details about git submodules see information :ref:`here<https://git-scm.com/book/en/v2/Git-Tools-Submodules>`.
   

.. ################################################
.. _install-BacterialTyper-source:

Install BacterialTyper from source
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Install Dependencies
""""""""""""""""""""

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
