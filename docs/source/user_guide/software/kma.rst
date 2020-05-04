.. ########################
.. _KMA-description:
.. ########################

KMA
====

The k-mer alignment (KMA) software creates an alignment method that 
allows for direct alignment of raw reads against entire databases, 
without the need of similarity reduction. In order to facilitate this, 
KMA uses an extra mapping step where the template of each input sequence 
is found and scored with the ConClave algorithm.