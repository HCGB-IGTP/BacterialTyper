.. _installing:

Installation
============

.. note::

    If you wish to contribute to the project, it's recommended you
    :ref:`install the latest development version<install_from_source>`.


This is an installation guide for BacterialTyper. 

.. _virual-env-BacterialTyper:

Python environment
------------------

It is also recommended to install BacterialTyper and all its dependencies and modules under a python environment.

Create environment
^^^^^^^^^^^^^^^^^^
To do so, run :file:`config/main/python_environment-commands.sh`, which is shown below:

.. include:: ../../../config/main/python_environment-commands.sh
   :literal:
   
Execute the file from the directory :file:`BacterialTyper/` as::

   sh config/main/python_environment-commands.sh
   
Activate environemt
^^^^^^^^^^^^^^^^^^^

Before executing any BacterialTyper module or script, it is convenient to activate the environment.

To do so, run :file:`config/main/activate_environment.sh`, which is shown below:

.. include:: ../../../config/main/activate_environment.sh
   :literal:

Execute the file from the directory :file:`BacterialTyper/` as::
   sh config/main/activate_environment.sh



.. _install_from_pip

Installing from pip
-------------------

Type::
   python -m pip install BacterialTyper


.. _install_from_source:

Installing from source
----------------------

If you are interested in contributing to BacterialTyper development, running the latest source code, 
or just like to build everything yourself, it is not difficult to build it from source. 

Grab the latest tar.gz release file from the PyPI files page, or if you want to develop BacterialTyper or just 
need the latest bugfixed version, grab the latest git version, and see Install from source.

Type::
   git clone https://github.com/JFsanchezherrero/BacterialTyper.git 

There are several extra python module requirements that are needed. 
They are listed in :file:`config/main/python_requirements.txt`, which is shown below:

.. include:: ../../../config/main/python_requirements.txt
   :literal:
   
Additionally, several third-party software packages are also required. 
They are listed in :file:`config/main/DEPENDENCIES`, which is shown below:

.. include:: ../../../config/main/DEPENDENCIES
   :literal:
   
