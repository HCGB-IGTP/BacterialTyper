#!/usr/bin/env python3
##########################################################
## Jose F. Sanchez                                      ##
## Copyright (C) 2019 Lauro Sumoy Lab, IGTP, Spain      ##
##                                                      ## 
##########################################################
##                                                      ## 
## Original version credit for:                         ## 
## https://gist.github.com/aziele/e060d1c0d1ef2acb4bd67e2fc0f165a6
## https://www.biostars.org/p/253984/
## 
## New implementations: add qlen and slen and           ## 
##                 control for aln & length thresh.     ## 
##                                                      ## 
##########################################################

"""
Script that parses BLAST output lists (in -outfmt '6 std qlen slen' format).

The original version was retrieved from:  

- https://gist.github.com/aziele/e060d1c0d1ef2acb4bd67e2fc0f165a6

- https://www.biostars.org/p/253984/

New implementation required query length (qlen) and subject lenght (slen) in the BLAST format input file. It allows to parse for alignment & length threshold. 

"""

from itertools import groupby

class Hsp:
    """Store information about single HSP in an alignment hit. 
    
    :param qid: Query Id
    :param sid: Subject Id
    :param pident: Percentage of identical matches
    :param length: Alignment length
    :param mismatch: Number of mismatches
    :param gaps: Total number of gaps
    :param qstart: Start of alignment in query
    :param qend: End of alignment in query
    :param sstart: Start of alignment in subject
    :param send: End of alignment in subject
    :param evalue: Expect value
    :param bitscore: Bit score
	:param qlen: Query length
	:param slen: Subject length
    """

    def __init__(self,entry):
        bt_fields = entry.split('\t')
        self.qid = bt_fields[0]
        self.sid = bt_fields[1]
        self.pident = float(bt_fields[2])
        self.length = int(bt_fields[3])
        self.mismatch = int(bt_fields[4])
        self.gaps = int(bt_fields[5])
        self.qstart = int(bt_fields[6])
        self.qend = int(bt_fields[7])
        self.sstart = int(bt_fields[8])
        self.send = int(bt_fields[9])
        self.evalue = float(bt_fields[10])
        self.bitscore = float(bt_fields[11])
        self.qlen = float(bt_fields[12])
        self.slen = float(bt_fields[13])

    def _format(self):
        l = [self.qid, 
             self.sid, 
             '{0:.2f}'.format(self.pident),
             '{0}'.format(self.length),
             '{0}'.format(self.mismatch),
             '{0}'.format(self.gaps),
             '{0}'.format(self.qstart),
             '{0}'.format(self.qend),
             '{0}'.format(self.sstart),
             '{0}'.format(self.send),
             '{:.1E}'.format(self.evalue),
             '{0}'.format(self.bitscore),
             '{0}'.format(self.qlen),
             '{0}'.format(self.slen)
            ]
        return l

    def format(self):
        return '\t'.join(self._format())

    def __str__(self):
        f = "{0}\t{1}".format(self.qid, self.sid)
        return f

class BlastRecord:
    """Object representing a Blast Record. 
    
    :param qid: Query sequence id
    :param hits: Blast hits
    """

    def __init__(self, qid=None, hits=None):
        """Initialize Blast Record instance"""
        self.qid = qid
        self.hits = hits

    def evalue_cutoff(self, evalue):
        """Filter HSPs by given e-value."""
        l = []
        for hit in self.hits:
            hsps = []
            for hsp in hit:
                if hsp.evalue <= evalue:
                    hsps.append(hsp)
            if hsps: l.append(hsps)
        self.hits = l

    def best_hits(self):
        """Return list of first hits that obtain the 
           same score.
        """
        l = []
        if not self.hits: return l
        else:
            max_score = self.hits[0][0].bitscore
            for hit in self.hits:
                if hit[0].bitscore >= max_score:
                    l.append(hit)
                else:
                    break
            return l

    def best_hsps(self):
        """Return list of first HSPs that obtain the 
           same score.
        """
        l = []
        if not self.hits: return l
        else:
            max_score = self.hits[0][0].bitscore
            for hit in self.hits:
                if hit[0].bitscore >= max_score:
                    l.append(hit[0])
                else:
                    break
            return l

    def best_hsps_except_query(self):
        """Return list of first HSPs that obtain the 
           same score.
        """
        l = []
        if not self.hits: return l
        else:
            max_score = None
            for hit in self.hits:
                if hit[0].sid != self.qid:
                    if not max_score or hit[0].bitscore>=max_score:
                        max_score = hit[0].bitscore
                        l.append(hit[0])
                    else:
                        break
            return l

    def _format(self):
        l = []
        for hit in self.hits:
            for hsp in hit:
                l.append(hsp._format())
        return l        

    def format(self):
        """Return output string of BLAST record"""
        l = []
        for hit in self.hits:
            for hsp in hit:
                l.append(hsp.format())
            l.append('')
        return "\n".join(l).strip()

    def __str__(self):
        """Return output string of BLAST record"""
        return self.format()


#This is a generator function!
def parse(handle, eval_thresh=10, aln_thresh=0, length_thresh=0):
     """Generator function to iterate over Blast records.
     
     :param handle: input file handle containg Blast tabular output -outfmt '6 std qlen slen'
     :param eval_thresh: E-value cutoff for Blast results. (1-eX)
     :param aln_thresh: Alingment cutoff (%)
     :param length_thresh: Cutoff for a given length provided (bp)
     
     :type handle: string
     :type eval_thresh: integer
     :type aln_thresh: integer
     :type length_thresh: integer
     
     :returns: BlastRecord objetct.
     
     """
     for qid, blasts in groupby(handle, lambda l: l.split()[0]):
         hits = []
         prev_sid = False
         for sid, hsps in groupby(blasts, lambda l: l.split()[1]):
             hsps_temp = []
             for line in hsps:
#                 line = line.functions.decode("utf-8")
                 hsp = Hsp(line)
                 
                 if (hsp.qid == hsp.sid):
                 	continue ## discard reporting autohit
                 
                 aln = (int(hsp.length)/int(hsp.qlen))*100
                 if aln >= aln_thresh and hsp.evalue <= eval_thresh and hsp.qlen > length_thresh:
                     hsps_temp.append(hsp)
             if hsps_temp: hits.append(hsps_temp)
         yield BlastRecord(qid=qid, hits=hits)



