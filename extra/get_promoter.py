from Bio import SeqIO
import sys

## read file containing proteins
if len(sys.argv)!=5:
    print()
    print ("## Usage protein_translator.py")
    print ("python %s gbf_file gene basePairs sampleName\n" %sys.argv[0])
    exit()

## arguments
gbf_file = sys.argv[1]
geneOfInterest = sys.argv[2]
basePairs = sys.argv[3]
sampleName = sys.argv[4]

for rec in SeqIO.parse(gbf_file, "genbank"):
	ID = rec.id
	SEQ = rec.seq

	## loop through features
	for feature in rec.features:
		if feature.type=="gene":
			qualif = feature.qualifiers
			for keys, values in qualif.items():
				#print (keys)
				#print (values)
				if values[0]==geneOfInterest:
					#print (feature)
					#print (ID)
					#print (feature.strand)
					
					if int(feature.strand) > 0:
						#print ("Start promoter: " + str(feature.location.nofuzzy_start-int(basePairs)))
						#print ("Start gene: " + str(feature.location.nofuzzy_start))
						#print ("End gene: " + str(feature.location.nofuzzy_end))
						promoter_seq = SEQ[feature.location.nofuzzy_start-int(basePairs):feature.location.nofuzzy_start]
						#, feature.location.nofuzzy_end           
					else:
						#print ("Start promoter: " + str(feature.location.nofuzzy_end+int(basePairs)))
						#print ("Start gene: " + str(feature.location.nofuzzy_end))
						#print ("End gene: " + str(feature.location.nofuzzy_start))
						promoter_seq = SEQ[feature.location.nofuzzy_end : feature.location.nofuzzy_end +int(basePairs)].reverse_complement()

					## print seq
					print (">" + sampleName + " promoter_" + basePairs + "_" + geneOfInterest)
					print (promoter_seq)  

