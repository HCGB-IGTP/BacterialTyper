************
Main modules
************

.. only:: html

    :Version: |version|
    :Date: |today|

There are multiple modules available that perform several :ref:`steps<pipeline-scheme>` and generate multiple results. 

``BacterialTyper`` modules require several command-line arguments and options to run. There are a number of shared 
:ref:`arguments<shared-arguments>` among all modules and some others specific of each module and specified 
within each one.

.. toctree::
   :maxdepth: 1
   
   annot.rst
   assemble.rst
   citation.rst
   cluster.rst
   config.rst
   database.rst
   help_info.rst
   ident.rst
   MGE.rst
   metadata.rst
   prep.rst
   profile.rst
   phylo.rst
   qc.rst
   report_generation.rst
   run.rst
   test.rst
   trimm.rst
   version.rst
   
.. _shared-arguments:

Command-line shared arguments
*****************************

Here we include a brief description of the shared command-line arguments for some of ``BacterialTyper`` modules.
   
   **Mode:**
   
   --project         Project mode. Requires as ``--input`` a folder containing an initialized ``BacterialTyper`` project [Default].
   
   --detached        Isolated mode. ``--input`` is a folder containining fastq reads. Provide a unique path o several using ``--batch`` option
   
   **Input/Output:**
   
   --input  string         Folder containing a project or reads, according to the mode selected. Files could be ``.fastq/.fq`` or ``fastq.gz/.fq.gz``. See ``--help_format`` for additional details.
   
   --single_end         Single end files [Default OFF]. Default mode is paired-end.
   
   --batch        Provide this option if input is a file containing multiple paths instead a path.
   
   --in_sample string         File containing a list of samples to include (one per line) from input folder(s) [Default OFF].
   
   --ex_sample string         File containing a list of samples to exclude (one per line) from input folder(s) [Default OFF].
   
   **Options:**
   
   --threads int        Number of CPUs to use [Default: 2]
   
   **Additional information:**
   
   --debug        Show additional message for debugging purposes.


.. _project-organization:

Details of the BacterialTyper project folder
**************************************************

TODO