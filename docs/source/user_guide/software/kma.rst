.. ########################
.. _KMA-description:
.. ########################

KMA
====

The k-mer alignment (KMA) software :cite:`Clausen2018` creates an alignment method that 
allows for direct alignment of raw reads against entire databases, 
without the need of similarity reduction. KMA uses an extra mapping step where the 
template of each input sequence is found and scored with the ConClave algorithm.

KMA - Output Guide

Explanation of the columns

Template: shows the name of the template sequences

Score: is the global alignment score of the template

Expected: is the expected alignment score if all mapping reads where smeared over all templates in the database

Template length: is the template length in nucleotides

template_id is the percent identity of the found template, over the full template length.

template_coverage is percent of the template that is covered by the query.

query_id is the percent identity between the query and template sequence, over the length of the matching query sequence.

query_coverage is the length of the matching query sequence divided by the template length.

Depth: is the number of times the template has been covered by the query.

q_value: is the quantile from McNemars test, to test whether the current template is a significant hit.

p_value: is p-value corresponding to the obtained q_value.

See additional details here: https://cge.cbs.dtu.dk/services/KMA/output.php