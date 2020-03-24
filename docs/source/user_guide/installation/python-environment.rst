.. ############################
.. _virtual-env-BacterialTyper:
.. ############################

Python environment
==================

It is always recommended to install ``BacterialTyper`` and all its dependencies and modules 
under a python environment. 


.. ###########
.. create-env:
.. ###########

Create environment
------------------
To do so, run :file:`env/create_python_environment.sh`, which is shown below:

.. include:: ../../../../env/create_python_environment.sh
   :literal:
   
Execute the file from the directory :file:`BacterialTyper/` as:

.. code-block:: sh

   sh env/create_python_environment.sh
   
.. ###########
.. _activate-env:
.. ###########

Activate environemt
-------------------

Before executing any ``BacterialTyper`` module or script, it is convenient 
to activate the environment.

To do so, run :file:`env/activate_environment.sh`, which is shown below:

.. include:: ../../../../env/activate_environment.sh
   :literal:

Execute the file from the directory :file:`BacterialTyper/` as:

.. code-block:: sh
  
   sh env/activate_environment.sh

Install modules
---------------

Once the ``BacterialTyper`` environment is created, install using ``pip`` ``BacterialTyper`` or any 
additional modules. 

*e.g.*: Install ``BacterialTyper``

.. code-block:: sh

   python -m pip install BacterialTyper

*e.g.*: install a list of requirements   

.. code-block:: sh
   
   pip install -r ./BacterialTyper/config/python_requirements.txt
   
.. attention:: Installing ``BacterialTyper`` using ``pip`` installs 
   all the python requirements. For additional third-party software and
   software required for ``BacterialTyper`` to work you will need to run the 
   ``config`` module.
