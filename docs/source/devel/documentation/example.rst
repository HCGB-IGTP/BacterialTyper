.. _example-doc:

.. function:: hello_world
   :noindex:

   Rename fasta sequences provided in file :file:`fasta_file` using id :file:`name`. 
   Save results in file :file:`new_fasta` provided.
   
   Check for id character lenght do not exceed 37 characters as it might be a limitiation 
   in further annotation and subsequent analysis. Read Prokka_ issue for further 
   details: https://github.com/tseemann/prokka/issues/337.

   :param fasta_file: Absolute path to fasta file.
   :param name: String to add every fasta sequence header.
   :param new_fasta: Name for the new fasta file (Absolute path).
   
   :type name: string
   :type fasta_file: string
   :type new_fasta: string
   
   :return: Path to tabular delimited file containing conversion from all to new id for each sequence.
   :warnings: Returns FAIL if name is >37 characters.
   
.. include:: ../../links.inc
   