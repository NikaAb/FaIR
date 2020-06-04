"""
Junction extraction with variable e !To Do! : integrate the parameters in order
to chooseing the e at the moment of execution.

"""

import sys
import difflib
from optparse import OptionParser
import numpy as np
import matplotlib.pyplot as plt

# =============================================================================
#               	
# =============================================================================	

def read_file(filename):
	filetoread = open(filename,"r")
	lines = filetoread.readlines()
	filetoread.close()
	return lines

def creat_fasta_dico(lines):
	seq = ""
	nom = ""
	fasta_dico = {}
	for l in lines:
		if l[0] == '>':
			if seq != "":
				fasta_dico[nom] = seq
			nom = l[1:-1]
			seq = ""
		else:
			seq = seq + l[:-1]
	if seq != "":
		fasta_dico[nom] = seq
	return fasta_dico

def find_position(lines):
	dico_position = {}
	for l in range(0,len(lines),4):
		if lines[l][0] == ">" :
			split = lines[l].split("\t")[1].split(" ")
			dico_position[lines[l].split("+")[0][1:-1]] = [split[1],split[2]]
	#print dico_position
	return dico_position

def extract_junction(fasta_dico,dico_position):
	dico_junction = {}
	for key in dico_position.keys():
		start = int(dico_position[key][0])
		stop = int(dico_position[key][1])
		dico_junction[key] = fasta_dico[key][start :stop].lower()
	return dico_junction

def write_junction_to_fasta(dico_junction,output_name) :
	filetowrite = open(output_name,"w")
	for key in dico_junction.keys():
		seqId = ">" + str(key) + "\n"
		filetowrite.write(seqId)
		seq = str(dico_junction[key]) + "\n" 
		filetowrite.write(seq)
	filetowrite.close()
	return 0
#=============================================================================#

def main():
    usage = "python unction_extraction.py -v <vidjil output> -f <fasta file> -o <Output_Fasta_file>\n "
    parser = OptionParser(usage)
    parser.add_option("-v", "--vidjil_output", dest="vidjil_output",
          help="the yaml file containing the partitions")
    parser.add_option("-f", "--Fasta_file", dest="Fasta_file",
          help="the name for the Fasta file to read")
    parser.add_option("-o", "--Output_Fasta_file", dest="Output_Fasta_file",
          help="the name for the Fasta file to write")
    (options, args) = parser.parse_args()
    if len(sys.argv) != 7:
    	parser.error("incorrect number of arguments")
    vidjil_file = options.vidjil_output
    Fasta_file = options.Fasta_file
    Output_Fasta_file = options.Output_Fasta_file
    vidjil = read_file(vidjil_file)

    fasta = read_file(Fasta_file)
    fasta_dico = creat_fasta_dico(fasta)
    dico_position = find_position(vidjil)
    dico_junction = extract_junction(fasta_dico,dico_position)
    write_junction_to_fasta(dico_junction,Output_Fasta_file)


#=============================================================================#

if __name__ == "__main__":
    main()
