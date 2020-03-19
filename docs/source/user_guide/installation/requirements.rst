.. ########################
.. _python-modules-required:
.. ########################

Python modules
--------------

There are several extra python module requirements that 
are needed and are summarized in the following table. 

.. csv-table::
   :header: "Package", "Version"
   :file: ../../data/python_requirements_summary.csv

These modules might have extra dependencies. The list of 
all modules required are listed in :file:`config/main/python_requirements.txt`.

.. include:: ../../../../config/main/python_requirements.txt
   :literal:

Although most of these dependencies will be fulfilled during the 
installation with ``pip``, you might be interested in installing them yourself. 
To install these dependencies see the appropiate section 
within :ref:`install from source<install-bacterialtyper-source>`, section 
:ref:`Install python modules<install-python-modules>`.

.. ########################
.. _perl_packages:
.. ########################

Perl packages
-------------

Although we tried to focus on Python_ as the main programming 
language for BacterialTyper some third-party software rely on perl 
and so some perl dependencies are required. 

They are listed in file 
:file:`config/main/perl_lib_dependencies.csv`, which is shown below:

.. csv-table::
   :header: "Package", "Version"
   :file: ../../data/perl_lib_dependencies_summary.csv

To install these dependencies see the appropiate section 
within :ref:`install from source<install-bacterialtyper-source>`, section 
:ref:`Install perl packages<install-perl-packages>`.

.. ######################
.. _soft-dependencies:
.. ######################

Software dependencies
---------------------

Also, several software packages are also required. They are listed in
:file:`config/main/DEPENDENCIES`, which is shown below:

.. csv-table::
   :header: "Third Party software"
   :file: ../../../../config/main/DEPENDENCIES

For additional details on how to install these dependencies see the 
appropiate section within :ref:`install from source<install-bacterialtyper-source>`, 
section :ref:`Install Software Dependencies<install-soft-deps>`.

.. ######################
.. _third-party_software:
.. ######################

Third party software
--------------------

Within the BacterialTyper project there is an additional module named 
IslandPath-DIMOB (IslandPath_ :cite:`Bertelli2017`) that we distribute 
along BacterialTyper. See additional details in the API section for 
:ref:`Third-party Software<third-party-soft>`.



.. #### Include links
.. include:: ../../links.inc