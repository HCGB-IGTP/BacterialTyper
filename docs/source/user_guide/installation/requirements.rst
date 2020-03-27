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

.. ########################
.. _perl_packages:
.. ########################

Perl packages
-------------

Although we tried to focus on Python_ as the main programming 
language for ``BacterialTyper`` some third-party software rely on perl 
and so some perl dependencies are required. 

They are listed in file 
:file:`BacterialTyper/config/perl/perl_lib_dependencies.csv`, which is shown below:

.. csv-table::
   :header: "Package", "Version"
   :file: ../../data/perl_lib_dependencies_summary.csv

To install these dependencies there is shell script available for the 
perl package modules installation (file :file:`BacterialTyper/config/perl_lib_installer.sh`). 

.. include:: ../../../../BacterialTyper/config/perl/perl_lib_installer.sh
   :literal:
 
Execute the file from the main directory as:

.. code-block:: sh

   sh BacterialTyper/config/perl/perl_lib_installer.sh

.. ######################
.. _soft-dependencies:
.. ######################

Software dependencies
---------------------

Also, several software packages are also required. They are listed in
:file:`BacterialTyper/config/dependencies.csv`, which is shown below:

.. csv-table::
   :header-rows: 1 
   :file: ../../../../BacterialTyper/config/dependencies.csv

Most of the software are common software that any person doing bioinformatics should have, so
you might have already available some of the software.

To test if any dependency is missing, you can execute the script: :file:`BacterialTyper/config/check-dependencies.sh`

.. code-block:: sh

   sh BacterialTyper/config/check-dependencies.sh

Once you identified the missing dependencies and minimum versions required you can either install them and 
set them available within your ``$PATH`` or you can execute the python script :file:`BacterialTyper/config/install_dependencies.py`.

.. ######################
.. _third-party_software:
.. ######################

Third party software
--------------------

Within the ``BacterialTyper`` project there is an additional module named 
IslandPath-DIMOB (IslandPath_ :cite:`Bertelli2017`) that we distribute 
along ``BacterialTyper``. See additional details in the API section for 
:ref:`Third-party Software<third-party-soft>`.

.. #### Include links
.. include:: ../../links.inc