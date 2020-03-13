.. ################################################
.. _faqs:


Frequently Asked Questions (FAQs)
*********************************

.. contents::

This is a collection of FAQs for BacterialTyper tutorial and interpretation of results. 

.. note::

    Please read further information on each section following the links provided or
    the main :ref:`User Guide for BacterialTyper<users-guide-index>`.


- Why read length values vary for each sample?

This is generally a result of trimming by default in the raw data output of the MiSeq run. The software bcl2fastq :cite:`bcl2fastq`
performs basecalling and writes the sequence and quality scores into two FASTQ files per sample (one for R1 and one
for R2), after separating demultiplexing mixed samples from the library pool by bar code indices. It also performs the
trimming of adapter sequence in the reads beyond the genomic insert. If so, the corresponding base calls beyond the
match will be changed to N in the resultant FASTQ file or trimmed, generating variable lengths for each read.

- Why is there a drop in quality with read length?

As seen in figure 2, there is a progressive drop in quality with read length, less pronounced in R1 than in R2, causing R2 reads with lower mean base
quality :cite:`Tan2019`, which is expected from 300 nt long paired end reads


.. bibliography:: ../bib/references.bib
.. include:: ../links.inc