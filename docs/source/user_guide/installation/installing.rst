.. ########################
.. _installing:
.. ########################

Installation
************

.. toctree::
   :hidden:
   
   requirements.rst
   python-environment.rst
      
.. contents::

This is an installation guide for the ``BacterialTyper`` pipeline. 

``BacterialTyper`` is a pipeline composed of multiple lines of code 
that calls and automatizes the analysis from multiple bioinformatic 
tools. Several dependencies (python, perl, software dependencies 
and third-party software) are required for the ``BacterialTyper`` analysis. 

You first need to get ``BacterialTyper`` (from different sources available): 
you can either install ``BacterialTyper`` latest version using `Python pip`_ 
or if you want to contribute or see additional details of the project, it's 
recommended you :ref:`install the latest development version<install-from-source>`.

Once you get the code, before running ``BacterialTyper`` you must make sure 
you have all the dependencies fulfilled from section 
:ref:`Requirements and dependencies<Requirements-dependencies>` either using the 
``config`` module or different scripts available for each type of requirement.


.. ########################
.. _Requirements-dependencies:
.. ########################

Requirements and dependencies
=============================

We include details on the different modules required here:

.. toctree:: 
   requirements.rst

.. ########################
.. _install-BacterialTyper:
.. ########################

Installing BacterialTyper
=========================

If you want to run a stable version of ``BacterialTyper``, install it using 
``pip`` package following instructions in :ref:`Installing from pip<install-from-pip>`.
On the other hand if you are interested in contributing to ``BacterialTyper`` development, 
running the latest source code, or just like to build everything yourself, you
should follow the :ref:`Installing from source<install-from-source>` instructions. 

Additionally, there are a number of dependencies that might be necessary to 
install or check within your system in both cases. Choose the appopiate choice 
according to your intalling ``BacterialTyper`` option selected.

  
.. ##################
.. _install-from-pip:
.. ##################

Installing from pip
-------------------

The fastest and most straightforward solution to install the latest 
release (or any desired version) of ``BacterialTyper`` is to install it (and 
all its python dependencies), using `Python pip`_. 

First, we encourage you to create a virtual environment before installing 
``BacterialTyper``. See  section for :ref:`Python environment<virtual-env-BacterialTyper>` 
for further details.

To install ``BacterialTyper`` using ``pip`` type:

.. code-block:: sh

   python3 -m pip install BacterialTyper
   
Follow additional `pip installing instructions`_ to learn about installing packages.

Also, as ``BacterialTyper`` relies in multiple dependencies, external perl packages, 
software and third-party software, we encourage you to once you install ``BacterialTyper`` 
using ``pip`` you check for dependencies using the ``BacterialTyper config`` module. See
details in ``BacterialTyper config`` module :ref:`section<config-description>`.
 
.. ##################
.. _install-from-source:
.. ##################

Installing from source
----------------------

Under some circumstancies (develop, bugs fixed, etc) you might be interested 
in obtaining the latest code version. Take into account, that you will need to install 
dependencies and fulfill requirements to have a working distribution. 

.. ##############
.. _get-git-code:
.. ##############

Get source code
^^^^^^^^^^^^^^^

The ``BacterialTyper`` project uses git_ as a version control system. To get the code, 
you can grab the latest version from the `BacterialTyper github`_ website, and 
follow the :ref:`Install BacterialTyper from source<install-bacterialtyper-source>` instructions.

Using the command-line, check you have a working distribution of ``git`` by typing 
``git --help`` or install it by typing:

.. code-block:: sh

   sudo apt update
   sudo apt upgrade
   sudo apt install git

Once you have ``git`` installed and working change directory to a folder of interest 
where the ``BacterialTyper`` project would be download. Type:

.. code-block:: sh

   git clone https://github.com/JFsanchezherrero/BacterialTyper.git 

Additionally, ``BacterialTyper`` contains a :ref:`third-party<third-party-soft>` project, 
IslandPath-DIMOB (IslandPath_ :cite:`Bertelli2017`) that ``BacterialTyper`` developers
forked, updated and adapted. See additional details on this project :ref:`here<islandPath-DIMOB-forked>`.

Then, if you are downloading ``BacterialTyper`` code from source, you need to retrieve 
the submodule latest version and install it accordingly. To do so, type:

.. code-block:: sh
   
   git clone --recurse-submodules https://github.com/JFsanchezherrero/BacterialTyper.git

On the other hand, you can also do it manually. In the directory containing ``git`` 
``BacterialTyper`` project you would find a ``third_party`` folder. Change directory to 
there and type:

.. code-block:: sh

   cd BacterialTyper/third_party
   git submodule init
   git submodule update
   
For additional details about ``git`` submodules see information here: 
https://git-scm.com/book/en/v2/Git-Tools-Submodules

.. ###############################
.. _install-bacterialtyper-source:
.. ###############################

Install BacterialTyper from source
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Once you have ``BacterialTyper`` source code available you only need to include
the ``BacterialTyper`` folder and main script within your path. 

.. code-block:: sh

   export PYTHONPATH=$PYTHONPATH":"$PWD"/BacterialTyper"

   export PATH=$PATH":"$PWD"/BacterialTyper.py"

Take into account that before running ``BacterialTyper`` you have to make sure you have all the 
dependencies fulfilled from section :ref:`Requirements and dependencies<Requirements-dependencies>`.
You can either install them yourself, use appropiate scripts for this purpose or use the ``BacterialTyper config``
module to check, update and install all dependencies required.


.. include:: python-environment.rst

.. ###########
.. include:: ../../links.inc
