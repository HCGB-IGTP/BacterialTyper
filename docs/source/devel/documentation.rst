.. _documentation:

Documentation
=============

Documentation is written using Sphinx_, a python documentation system built using 
reStructuredText_ (ReST; ``.rst``). The docs configuration contains both 
ReST files that contain pages in the documentation and configuration files for Sphinx_.

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


Build the docs
--------------

Instructions to build the documentation for developer purposes.

All documentation is built from the :file:`BacterialTyper/docs/` directory. 

We will follow the rules for the documentation generated for the Matplotlib documentation configuration_.
 
.. ###############################################
.. _docs-organization-folder:

.. ###############################################
.. _installing-dep-build-docs:

Installing dependencies
^^^^^^^^^^^^^^^^^^^^^^^

You may need to


.. ###############################################
.. _building-docs-guide:

Building documentation
^^^^^^^^^^^^^^^^^^^^^^

To build the docs...


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

Sphinx_ allows internal references_ between documents.

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


.. ###############################################
.. References
.. ###############################################
.. _Sphinx: http://www.sphinx-doc.org/en/master/
.. _reStructuredText: http://docutils.sourceforge.net/rst.html
.. _configuration: https://matplotlib.org/devel/documenting_mpl.html
.. _api: https://en.wikipedia.org/wiki/Application_programming_interface
.. _references: https://www.sphinx-doc.org/en/stable/usage/restructuredtext/roles.html
