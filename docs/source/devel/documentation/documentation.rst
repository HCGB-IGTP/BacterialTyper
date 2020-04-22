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

The documentation folder :file:`BacterialTyper/docs/source` organization is distributed 
in multifple files, folders and subfolders:

* :file:`api` - Placeholders to automatically generate the application programming interface (API_) documentation.

* :file:`devel` - Developers guidelines to contribute to the project.

* :file:`faq` - Frequently asked questions (FAQs).

* :file:`glossary` - Glosary entries and abreviations.

* :file:`images` - Images to include in the documentation.

* :file:`tutorial` - Tutorial examples for the ``BacterialTyper`` modules.

* :file:`user_guide` - The user guide and documentation to interpret results and analysis.

* :file:`index.rst` - The top level index document for the documentation.

* :file:`conf.py` - Sphinx configuration parameters.

* :file:`Makefile` and :file:`make.bat` - Entry points for building the docs.

* :file:`_static` - Used by the sphinx build system.

* :file:`_templates` - Used by the sphinx build system.

The ``.rst`` files are kept in :file:`user_guide`, :file:`devel`, 
:file:`api`, :file:`glossary` and :file:`tutorial`. 

The main entry point is :file:`index.rst`, which contains links
for users guide, developers guide, api reference, and FAQs. 
The documentation suite is built as a single 
document in order to make the most effective use of cross referencing.

There are also ``.rst`` files that are contained in :file:`api/modules` 
and :file:`api/scripts` that are automatically generated from the docstrings 
of the functions in ``BacterialTyper`` scripts and main modules. These sources consist 
of python scripts that have ReST documentation built into their comments. See section
:ref:`api-docstrings` for details. 

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

To build the docs, you will need to install several python modules because 
documentation is generated from reStructuredText_ (ReST) using the 
Sphinx_ documentation generation tool. 

.. attention::

  * You will need a minimal working ``LaTeX`` distribution.

There are several extra ``python`` requirements that are needed to build the documentation. 
They are listed in :file:`docs/config/doc-requirements.txt`, which is shown below:

.. include:: ../../../config/doc-requirements.txt
   :literal:

You will need a ``BacterialTyper`` (and dependencies) working distribution included in your
``$PYTHONPATH`` (activate your virtual environment as mentioned :ref:`here<activate-env>`) and then 
additionally install documentation requirements using ``pip``.

.. code-block:: sh

   pip install docs/config/doc-requirements.txt

Also, there is a ``perl`` module required for the documentation generation:

.. include:: ../../../config/perl_packages.txt
   :literal:

We encorage to install this module using CPAN (https://www.cpan.org/) or any other ``perl`` package
installer of your interest.

Finally, in order to generate the ``python`` graph dependency images, the tool ``graphviz`` and 
an specific ``python`` module is required.

To install ``graphviz`` type:

.. code-block:: sh

   sudo apt install graphviz
   
The ``python`` module ``pyan`` needs to be installed using ``pip`` from ``github``. 
See https://github.com/ttylec/pyan for additional details. Type:

.. code-block:: sh

   pip install git+https://github.com/ttylec/pyan

.. ###############################################
.. _building-docs-guide:

Building documentation
^^^^^^^^^^^^^^^^^^^^^^
The documentation sources are found in the :file:`docs/source` directory in the trunk.
The configuration file for Sphinx is :file:`docs/source/conf.py`. It controls which
directories Sphinx parses, how the docs are built, and how the extensions are
used. 

To build the documentation in html format, cd into :file:`docs/source/` and run:

.. code-block:: sh

   make html

To delete built files. It may help if you get errors about missing paths or broken links.

.. code-block:: sh

   make clean

To generate a pdf file of the documentation.
   
.. code-block:: sh

   make latexpdf

.. _writing-rest-pages:

Writing ReST pages
------------------

Most documentation is either in the docstring of individual classes and methods, 
in explicit ``.rst`` files, or in examples and tutorials.

All of these use the reStructuredText_ (ReST) syntax. Users should look at the ReST documentation
for a full description. But some hints and conventions uses are describing as a quickstart
for creating documentation.

.. ###############################################
.. _internal-section-refs:

Referring to other documents and sections
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Sphinx_ allows `internal references`_ between documents.

Documents can be linked with the ``:doc:`` directive:

.. code-block:: rst

   See the :doc:`../../../info/info_index`

   See the installation user guide :doc:`../../../user_guide/installation/installing`

will render as:

  See the :doc:`../../../info/info_index`
  
  See the installation user guide :doc:`../../../user_guide/installation/installing`
  
Sections can also be given reference names and be referenced using the ``:ref:`` directive.  

For instance, see here an example from the :doc:`/user_guide/modules/assemble` link:

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

.. ###############
.. _hyperlink-ref:

Include hyperlinks
^^^^^^^^^^^^^^^^^^

To include an hyperlink in a word just type:

.. code-block:: rst
   
   Example: `Google <https://www.google.com>`_.

and this will create the following: 

Example: `Google <https://www.google.com/>`_

.. ###############
.. _citations-ref:

Include citations
^^^^^^^^^^^^^^^^^

The bibliography is included ``BibTex`` format in file :file:`docs/source/bib/references.bib`. 
The initiation and style of the bibliography is determined by file :file:`docs/source/bib/bib_index.rst`
and the directive ``.. bibliography::``.

To include a citation you will use the role ``:cite:`` and the tag associated to each entry.

And example of citation might be:

.. code-block:: rst
   
   Whatever text you want to reference :cite:`Deplano2018`.

That will render as:

Whatever text you want to reference :cite:`Deplano2018`

.. ##########################
.. _writing-docstrings:

Writing docstrings
------------------

Docstrings should conform to the reStructuredText_ (ReST) syntax guide. 

.. #################
.. _example-docstrings:

Example docstring
^^^^^^^^^^^^^^^^^

An example docstring looks like:

.. literalinclude:: ./example.py
   :language: python
     

And it displays like:

.. toctree:: 
   :hidden:
     
   example.rst

.. include:: ./example.rst

For other docstrings example visit https://matplotlib.org/devel/documenting_mpl.html#example-docstring

Check an example for python documentation here: https://thomas-cokelaer.info/tutorials/sphinx/index.html


.. #################
.. _api-docstrings:

API Docstrings
^^^^^^^^^^^^^^

Most of the API documentation is written in docstrings. These are comment
blocks in source code that explain how the code works. For each module and 
script available we have to create an ``.rst`` file in :file:`docs/api` 
directory. Then, using the ``:automodule:`` directive and ``:members:`` role
we automatically include every docstring within the source code available.

Using a shell loop we create an ``.rst`` file for each module and script:

- modules (file :file:`docs/source/devel/documentation/modules_api-docstrings.sh`):

.. literalinclude:: ./modules_api-docstrings.sh
   :language: sh
  
- scripts (file :file:`docs/source/devel/documentation/scripts_api-docstrings.sh`):

.. literalinclude:: ./scripts_api-docstrings.sh
   :language: sh

Python Graph images
^^^^^^^^^^^^^^^^^^^

Additionally, for each script or module of interest, we generate a graph representation 
of the different functions and relationships with other modules. 

Images are stored in  :file:`docs/source/images/python_graph`. We also include automatically 
these images generated for the corresponding ``.rst`` file using the previous shell script.

To generate the ``python`` graph image, we would employ the ``python`` module ``pyan`` 
and the ``graphviz`` tool. See https://github.com/ttylec/pyan for additional details.


.. ## Include linksReferences
.. include:: ../../links.inc