import sys
from Bio import SeqIO, Seq
from Bio.SeqRecord import SeqRecord

## read file containing proteins
if len(sys.argv)!=3:
	print()
	print ("## Usage protein_translator.py")
	print ("python %s gene_CDS_file outfile\n" %sys.argv[0])
	exit()


## arguments
genes_file = sys.argv[1]
out_file = sys.argv[2]

with open (genes_file) as in_handle:
           prot_recs = SeqIO.to_dict(SeqIO.parse(in_handle, "fasta"))

trans_prot = {}
for prot in prot_recs.items():
	ID= prot[1].id
	prot[1].seq = prot[1].seq.translate(table=11, to_stop=True)
	#trans_prot[ID]=SEQ

	print ("+ Tranlate DNA seq for " + ID)


## save output
with open(out_file, "w") as handle:
    SeqIO.write(prot_recs.values(), handle, "fasta") 
