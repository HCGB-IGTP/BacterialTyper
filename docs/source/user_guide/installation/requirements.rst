.. ########################
.. _python-modules-required:
.. ########################

Python modules
--------------

There are several extra python module requirements that 
are needed and are summarized in the following table. 

.. csv-table::
   :header: "Module", "Version"
   :file: ../../data/installation/python_requirements_summary.csv

These modules might have extra dependencies. Details of the list of 
all modules required are listed in :file:`BacterialTyper/config/python/python_requirements.txt`. 
And accessible here:

.. toctree:: 
   python-requirements.rst

Although these dependencies will be fulfilled during the ``BacterialTyper``
installation with ``pip``, you might be interested in installing them yourself. 

Using ``pip`` we can install them all at a glance. 

.. code-block:: sh

   pip install -r ./BacterialTyper/config/python/python_requirements.txt

But again, following installation recommendations, we encourage you to create and install them 
within a virtual environment (See section: :ref:`Python environment<virtual-env-BacterialTyper>` 
section for details).

You can test the presence of these ``python`` modules using the ``BacterialTyper config`` module. 
Once you identified the missing dependencies and minimum versions required you can either install them and 
set them available within your ``PYTHONPATH`` or environment or you can execute the ``BacterialTyper config`` 
with ``install`` option.


.. ########################
.. _perl_packages:
.. ########################

Perl packages
-------------

Although we tried to focus on ``Python`` as the main programming 
language for ``BacterialTyper`` there are some scripts that are in ``perl``. 
This is due to simplicity as this functions were previously
written in other projects and we just re-used the code. As time goes by, we would 
like to reduce the ``perl`` code to a minimum and rewrite it in ``python``.

The ``perl`` packages required within these scripts are core modules of ``perl`` and should
be available within your installation. These modules are listed in file 
:file:`BacterialTyper/config/perl/perl_dependencies.csv`, which is shown below:

.. csv-table::
   :header-rows: 1
   :file: ../../../../BacterialTyper/config/perl/perl_dependencies.csv


You can test the presence of these ``perl`` modules using the ``BacterialTyper config`` module. 

.. ######################
.. _soft-dependencies:
.. ######################

Software dependencies
---------------------

Also, several software packages are also required. They are listed in
:file:`BacterialTyper/config/software/dependencies.csv`, which is shown below:

.. csv-table::
   :header-rows: 1 
   :file: ../../../../BacterialTyper/config/software/dependencies.csv

Most of the software are common software that any person doing bioinformatics should have, so
you might have already available within your system.

You might need to have installed some basic libraries: the default ``java`` run time environment (jre), C-compiler 
and zlib development files:

.. code-block:: sh
	
	sudo apt install default-jre
	sudo apt-get install libz-dev
	sudo apt-get install build-essential
	
You can test for any missing software dependencies using the ``BacterialTyper config`` module. Once you 
identified the missing dependencies and minimum versions required you can either install them and 
set them available within your ``$PATH`` or you can execute the ``BacterialTyper config`` 
with ``install`` option.


.. ##########################
.. _third-party-requirement:
.. ##########################

Third party software
--------------------

Additionally, within the ``BacterialTyper MGE`` module there is an optional analysis
that determines the pathogenic islands within the genome. The tool that we
employ for this purpose is named IslandPath-DIMOB (IslandPath_ :cite:`Bertelli2017`) 
and we distribute it along ``BacterialTyper``. See additional details in the API section for 
:ref:`Third-party Software<third-party-soft>`.

This third party software requires multiple ``perl`` modules and in turn, depend on many others,
making the installation and distribution quite difficult. As we established this
analysis as supplementary, so far, we would not take care of fulfilling the multiple
requirements of this tool. This packages are listed in here:

.. csv-table::
   :header: "Package", "Version", "Type_module"
   :file: ../../data/installation/perl_dependencies_IslandPath_summary.csv

We encorage to install these modules using CPAN (https://www.cpan.org/) or any other ``perl`` package
installer of your interest.

To test the presence of these additional ``perl`` modules, you can also use the ``BacterialTyper config`` module.



.. #### Include links
.. include:: ../../links.inc