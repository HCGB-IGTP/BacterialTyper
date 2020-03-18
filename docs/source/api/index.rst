.. title:: API guide

API Overview
============

.. toctree::
   :hidden:
   
.. only:: html

    :Version: |version|
    :Date: |today|

.. contents:: :local:

This pipeline is composed of multiple modules and scripts that 
are separated in four directories.

- ``BacterialTyper/modules``
- ``BacterialTyper/scripts``
- ``BacterialTyper/data``
- ``BacterialTyper/other_tools``
   
The main BacterialTyper script is within the folder ``main/`` that integrates and 
connects all available modules and analysis.

.. toctree::
   :maxdepth: 2

   BacterialTyper.rst
   modules/index.rst
   scripts/index.rst
   other_tools.rst
   data.rst