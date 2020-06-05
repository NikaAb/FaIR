import sys
import math
import resource
import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
from collections import Counter
from optparse import OptionParser
from sklearn.cluster import KMeans
from sklearn.cluster import AffinityPropagation



#####################################################################
def read_file (nomFi):
	f=open(nomFi,"r")
	lines=f.readlines()
	f.close()
	return lines
#####################################################################
#read the fasta file of input sequences
def readFastaMul(nomFi):
	lines = read_file (nomFi)
	seq=""
	nom=""
	lesSeq={}
	outliers ={}
	for l in lines:
		if l[0] == '>':
			if seq != "":
				if len (seq) <= 450 and len(seq) >= 350:
					lesSeq[nom] = seq
				else :
					outliers[nom] = seq
			nom=l[1:-1]
			seq=""
		else:
			seq=seq+l[:-1]
	if seq != "":
		if len (seq) <= 450 and len(seq) >= 350:
			lesSeq[nom.rstrip()] = seq.rstrip()
		else :
			outliers[nom.rstrip()] = seq.rstrip()
	print "number of outliers : ", len(outliers)
	return lesSeq
#####################################################################
def extract_junction(seq_dico,V_end,J_start):
	junction_dico = {}
	for key in seq_dico.keys():
		junction = seq_dico[key][V_end : -(J_start + 1)]
		if len(junction) < 15 :
			print key 
		else :
			junction_dico[key] = junction 
	return junction_dico


def write_junction_to_fasta(junction_dico,output_name) :
	filetowrite = open(output_name,"w")
	for key in junction_dico.keys():
		seqId = ">" + str(key) + "\n"
		filetowrite.write(seqId)
		seq = str(junction_dico[key]) + "\n" 
		filetowrite.write(seq)
	filetowrite.close()
	return 0
####################################################################
def main():
	usage = "usage: seq_length_distribution.py -f FastaFile"
	parser = OptionParser(usage)
	parser.add_option("-f", "--FastaFile", dest="FastaFile",
	      help="read data from FILENAME")
	(options, args) = parser.parse_args()
	if len(sys.argv) != 3:
		parser.error("incorrect number of arguments")
	
	FastaFile = options.FastaFile
	V_end = 280
	J_start = 27
	seq_dico = readFastaMul(FastaFile)
	junction_dico = extract_junction(seq_dico,V_end,J_start)
	output_name = str(FastaFile.split(".")[0]) + "_junction.txt"
	write_junction_to_fasta(junction_dico,output_name) 
#####################################################################
if __name__ == "__main__":
	main()
  