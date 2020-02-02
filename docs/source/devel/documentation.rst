.. _documentation:

Documentation
=============

.. contents:: Contents

Documentation is written using Sphinx_, a python documentation system built using 
reStructuredText_ (ReST; ``.rst``). The docs configuration contains both 
ReST files that contain pages in the documentation and configuration files for Sphinx_.

The Sphinx_ website also contains plenty of `documentation details concerning ReST`_
markup and working with Sphinx in general.

.. ###############################################
.. _docs-organization-folder:

Organization
------------

The documentation folder :file:`BacterialTyper/docs` organization is distributed in multifple files, folders and subfolders:

* :file:`api` - Placeholders to automatically generate the application programming interface (api_) documentation.

* :file:`devel` - Developers guidelines to contribute to the project.

* :file:`faq` - Frequently asked questions (FAQs).

* :file:`glossary` - Glosary entries and abreviations.

* :file:`images` - Images to include in the documentation.

* :file:`tutorial` - Tutorial examples for the BacterialTyper modules.

* :file:`user_guide` - The user guide and documentation to interpret results and analysis.

* :file:`index.rst` - The top level index document for the documentation.

* :file:`conf.py` - Sphinx configuration parameters.

* :file:`Makefile` and :file:`make.bat` - Entry points for building the docs.

* :file:`_static` - Used by the sphinx build system.

* :file:`_templates` - Used by the sphinx build system.

The ``.rst`` files are kept in :file:`user_guide`, :file:`devel`, 
:file:`api`, :file:`glossary` and :file:`tutorial`. 

The main entry point is :file:`index.rst`, which pulls in 
the :file:`index.rst` file for the users guide, developers guide, 
api reference, and FAQs. The documentation suite is built as a single 
document in order to make the most effective use of cross referencing.

.. note::

   There are also ``.rst`` files that are contained in :file:`api/modules` 
   and :file:`api/submodules` that are automatically generated from the docstrings 
   of the functions in BacterialTyper scripts and main modules. These sources consist 
   of python scripts that have ReST documentation built into their comments. 

.. _buld-the-docs:

Build the docs
--------------

Instructions to build the documentation for developer purposes.

All documentation is built from the :file:`BacterialTyper/docs/` directory. 

We will follow the rules for the documentation generated for the `Matplotlib documentation configuration`_.
 
.. ###############################################
.. _installing-dep-build-docs:

Installing dependencies
^^^^^^^^^^^^^^^^^^^^^^^

To build the docs, you will need to install several python modules as 
documentation is generated from reStructuredText_ (ReST) using the 
Sphinx_ documentation generation tool. 

There are several extra requirements that are needed to build the documentation. 
They are listed in :file:`config/docs/doc-requirements.txt`, which is shown below:

.. include:: ../../../config/docs/doc-requirements.txt
   :literal:

.. note::
  * You'll need a minimal working LaTeX distribution for many examples to run.


.. ###############################################
.. _building-docs-guide:

Building documentation
^^^^^^^^^^^^^^^^^^^^^^
The documentation sources are found in the :file:`docs/` directory in the trunk.
The configuration file for Sphinx is :file:`docs/conf.py`. It controls which
directories Sphinx parses, how the docs are built, and how the extensions are
used. 

To build the documentation in html format, cd into :file:`docs/` and run:

.. code-block:: sh

   make html

To delete built files. It may help if you get errors about missing paths or broken links.

.. code-block:: sh

   make clean

To generate a pdf file of the documentation.
   
.. code-block:: sh

   make latexpdf


.. note::

   The ``SPHINXOPTS`` variable is set to ``-W --keep-going`` by default to build
   the complete docs but exit with exit status 1 if there are warnings.  To unset
   it, use

   .. code-block:: sh
   
      make SPHINXOPTS= html
   
   You can use the ``O`` variable to set additional options:
   
   * ``make O=-j4 html`` runs a parallel build with 4 processes.
   * ``make O=-Dplot_formats=png:100 html`` saves figures in low resolution.
   * ``make O=-Dplot_gallery=0 html`` skips the gallery build.
   
   Multiple options can be combined using e.g. ``make O='-j4 -Dplot_gallery=0'
   html``.
   
   On Windows, options needs to be set as environment variables, e.g. ``set O=-W
   --keep-going -j4 & make html``.

.. _writing-rest-pages:

Writing ReST pages
------------------

Most documentation is either in the docstring of individual classes and methods, 
in explicit ``.rst`` files, or in examples and tutorials.

All of these use the reStructuredText_ (ReST) syntax. Users should look at the ReST documentation
for a full description. But some specific hints and conventions Matplotlib
uses are useful for creating documentation.


.. ###############################################
.. _internal-section-refs:

Referring to other documents and sections
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Sphinx_ allows `internal references`_ between documents.

Documents can be linked with the ``:doc:`` directive:

.. code-block:: rst

   See the :doc:`/info/history`

   See the installation user guide :doc:`/user_guide/installing`

will render as:

  See the :doc:`/info/history`

  See the installation user guide :doc:`/user_guide/installing`
  
Sections can also be given reference names.  For instance from the
:doc:`/user_guide/modules/assemble` link:

.. code-block:: rst

   .. _assembly-workflow:

   Workflow
   --------

   The assemble module contains several main functions (See :doc:`assemble module <../../api/modules/assemble>` for additional details.)
   Below we show the workflow of the assembly process. 

and refer to it using the standard reference syntax:

.. code-block:: rst

   See :ref:`assembly-workflow`

will give the following link: :ref:`assembly-workflow`

.. note::

   To maximize internal consistency in section labeling and references,
   use hyphen separated, descriptive labels for section references. Since 
   underscores are widely used by Sphinx itself, use hyphens to separate words.

.. _writing-docstrings:

Writing docstrings
==================

Most of the API documentation is written in docstrings. These are comment
blocks in source code that explain how the code works.

All new or edited docstrings should conform to the `numpydoc docstring guide`_.
Much of the ReST syntax discussed above (:ref:`writing-rest-pages`) can be
used for links and references.  These docstrings eventually populate files in 
:file:`docs/api` directory and form the reference documentation for the
library.


Example docstring
-----------------

An example docstring looks like:

.. code-block:: python

   """
   def rename_fasta_seqs(fasta_file, name, new_fasta):
   
      Rename fasta sequences provided in file :file:`fasta_file` using id :file:`name`. Save results in file :file:`new_fasta` provided.
      
      Check for id character lenght do not exceed 37 characters as it might be a limitiation in further annotation and subsequent analysis. Read Prokka_ issue for further details: https://github.com/tseemann/prokka/issues/337.
      
      :param fasta_file: Absolute path to fasta file.
      :type fasta_file: string
      :param name: String to add every fasta sequence header.
      :type name: string
      :param new_fasta: Name for the new fasta file (Absolute path).
      :type new_fasta: string
      :return: Path to tabular delimited file containing conversion from all to new id for each sequence.
      :warnings: Returns FAIL if name is >37 characters.
      
      .. include:: ../../links.inc
      """     
   
See some other example in https://matplotlib.org/devel/documenting_mpl.html#example-docstring

Check and example for python documentation here: https://thomas-cokelaer.info/tutorials/sphinx/index.html

.. ## Include linksReferences
.. include:: ../links.inc