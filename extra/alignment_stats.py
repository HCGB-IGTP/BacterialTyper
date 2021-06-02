from Bio import SeqIO
from Bio.Align import AlignInfo
from Bio import AlignIO
import sys

## read file containing proteins
if len(sys.argv)!=2:
	print()
	print ("## Usage protein_translator.py")
	print ("python %s alignment\n" %sys.argv[0])
	exit()


## arguments
aln_file = sys.argv[1]

align = AlignIO.read(aln_file, 'fasta') # Whatever format your mutliple sequence alignment is in

print ("#########################")
print ("align")
print (align)
print ()

print ("#########################")
print ("Detail positions")
poly_positions = []
for i in range(1,len(align[1,:])):
    list_entries = list(set(align[:,i]))
    if "-" in list_entries:
        list_entries.remove('-')

    chars = "".join(list_entries)
    print (str(i) + ": " + chars)

    if len(list_entries) > 1:
       poly_positions.append(i)

print ()

print ("#########################")
print ("Polymorphisms")
print (poly_positions)
for record in align: 
    for j in poly_positions:
        print (record.id + ":" + str(j) + ":" + record.seq[j])

exit()

print ("#########################")
summary_align = AlignInfo.SummaryInfo(align)
print ("summary_align")
print (summary_align)
print ()

print ("#########################")
consensus = summary_align.dumb_consensus()
print ("consensus")
print (consensus)
print ()

print ("#########################")
my_pssm = summary_align.pos_specific_score_matrix(consensus,chars_to_ignore = ['X'])
print ("my_pssm")
print(my_pssm)
print ()





