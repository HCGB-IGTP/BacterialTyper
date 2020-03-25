.. ############################
.. _virtual-env-BacterialTyper:
.. ############################

Python environment
==================

It is always recommended to install ``BacterialTyper`` and all its dependencies and modules 
under a `python environment`_. 

By executing the ``env/create_python_environment.sh`` script, it creates an environment, installs
all python requirements and activates it for later usage.

You only have to create the environment once and you can populate with additional modules,
see details in :ref:`Install modules<install-modules-env>`. 

Make sure to deactivate the environment after finishing the exection of ```BacterialTyper``
or activate it before running it, following the instructions in :ref:`Deactivate environment<deactivate-env>`
or :ref:`Activate environment<activate-env>` sections respectively. 

.. ###########
.. create-env:
.. ###########

Create environment
------------------
To do so, run :file:`env/create_python_environment.sh`, which is shown below:

.. include:: ../../../../env/create_python_environment.sh
   :literal:

You will need as a requisite to have installed in your system ``python3-dev`` and ``python3-venv``libraries. 

Install them by typing:

.. code-block:: sh

   sudo apt install python3-dev
   sudo apt install python3-venv

Execute the file from the main directory as:

.. code-block:: sh

   source env/create_python_environment.sh
   
Within this create environment script, there is a command to install all python
requirements for ``BacterialTyper``:

.. code-block:: sh

   pip install -r BacterialTyper/config/python_requirements.txt
    
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
  
   source env/activate_environment.sh
   

.. ###########
.. _deactivate-env:
.. ###########

Deactivate environemt
---------------------

After finished the execution of any ``BacterialTyper`` module or script, it is convenient 
to deactivate the environment. You can just close the terminal but it would be more appopiate
to conveniently deactivate the environment first.

To do so, execute the following command:

.. code-block:: sh
  
   deactivate



   
.. ###########
.. _install-modules-env:
.. ###########

Install modules
---------------

Once the ``BacterialTyper`` environment is created, install using ``pip`` ``BacterialTyper`` or any 
additional modules. 

*e.g.*: Install ``BacterialTyper``

.. code-block:: sh

   python -m pip install BacterialTyper

*e.g.*: install a list of requirements   

.. code-block:: sh
   
   pip install -r BacterialTyper/config/python_requirements.txt
   
.. attention:: Installing ``BacterialTyper`` using ``pip`` installs 
   all the python requirements. For additional third-party software and
   software required for ``BacterialTyper`` to work you will need to run the 
   ``config`` module.
