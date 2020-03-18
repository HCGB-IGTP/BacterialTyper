.. ########################
.. _python-modules-required:
.. ########################

Python modules
--------------

There are several extra python module requirements that are needed and are summarized in the following table. 

.. csv-table::
   :header: "Package", "Version"
   :file: ../../data/python_requirements_summary.csv

These modules might have extra dependencies. The list of all modules required are listed in :file:`config/main/python_requirements.txt`.

.. include:: ../../../../config/main/python_requirements.txt
   :literal:

Most of these dependencies will be fulfilled during the installation with ``pip``.

To install these dependencies see the appropiate section within :ref:`install-BacterialTyper`.

.. ########################
.. _perl_packages:
.. ########################

Perl packages
-------------

Although we tried to focus on Python_ as the main programming language for BacterialTyper some third-party software
rely on perl and so some perl depencies are required. They are listed in file :file:`config/main/perl_lib_dependencies.csv`, 
which is shown below:

.. csv-table::
   :header: "Package", "Version"
   :file: ../../data/perl_lib_dependencies_summary.csv

To install these dependencies see the appropiate section within :ref:`install-BacterialTyper`.

.. ######################
.. _third-party_software:
.. ######################

Third party software
--------------------

Also, several third-party software packages are also required. 
They are listed in :file:`config/main/DEPENDENCIES`, which is shown below:

.. csv-table::
   :header: "Third Party software"
   :file: ../../../../config/main/DEPENDENCIES

For additional details on how to install these dependencies see the appropiate section within :ref:`install-BacterialTyper`.


.. #### Include links
.. include:: ../../links.inc