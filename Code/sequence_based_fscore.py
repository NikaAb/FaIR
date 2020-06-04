"""
This script compares real clusters and predicted clusters content and measure 
precision and recall.

The comparison is based on method explained in method file. 

True cluster's file and predicted cluster's file have the same format:

1	a b
2	c e
3	d

The first column shows the cluster label and the second one represents sequences 
in that cluster, separated by space.
Two columns are seperated by \t.
"""


import sys
from collections import Counter 
from optparse import OptionParser	
#===================================================================================
def read_output_file(filename):
	f=open(filename,"r")
	lines=f.readlines()
	f.close()
	return lines
#===================================================================================
def creat_hashtrue(lines):
	HashTrue={}	
	label_seq_True = {}	
	for l in lines:
		label=l.split("\t")[0]
		sequences=l.split("\t")[1].split(" ")
		HashTrue[label]=len(sequences)
		for s in sequences :
			label_seq_True[s.rstrip()] = label
	return HashTrue, label_seq_True
#===================================================================================

def read_hashpredict(lines,label_seq_True):
	#print label_seq_True
	HashPredict={} # contains lables of each cluster
	seq_dico={}  # contains sequences of each cluster
	for l in lines: 
		label=l.split("\t")[0]

		sequences=l.split("\t")[1].rstrip().split(" ")
		seq_dico[label]=sequences
		for s in sequences:
			if label in HashPredict.keys():
				if s.rstrip() in label_seq_True.keys():
					HashPredict[label].append(label_seq_True[s.rstrip()])
				else :
					print "not in true_cluster : ", s.rstrip()
			else:
				if s.rstrip() in label_seq_True.keys():
					HashPredict[label]=[label_seq_True[s.rstrip()]]
				else :
					print "not in true_cluster : ", s.rstrip()
		#print HashPredict
		if label in HashPredict.keys():
			HashPredict[label]=dict(Counter(HashPredict[label]))
	return HashPredict,seq_dico

#===================================================================================
def counter(seq_dico,HashTrue,HashPredict,label_seq_True):
	TP,FN,FP=0,0,0
	for key in seq_dico.keys():
		tp,fn,fp=0,0,0
		for s in seq_dico[key]:
			if s.rstrip() in label_seq_True.keys():
				tp=int(HashPredict[key][label_seq_True[s.rstrip()]])-1
				fn=int(HashTrue[label_seq_True[s.rstrip()]])-1-tp
				for c in HashPredict[key].keys():
					if c !=label_seq_True[s.rstrip()]:
						fp=int(HashPredict[key][c])		
				TP+=tp
				FN+=fn
				FP+=fp
	return TP,FN,FP
#===================================================================================
def Classiq_mesure(TP,FN,FP):
	pre = TP/float(TP+FP)
	rec = TP/float(TP+FN)
	print "Recall : ", rec
	print "Precision : ", pre
	print "F1-score : ", round(2*pre*rec/(pre + rec),2)
	return 0
#===================================================================================
def precision_recall(Predicted_File,True_file):
	predicted_lines=read_output_file(Predicted_File)
	real_lines=read_output_file(True_file)
	HashTrue, label_seq_True = creat_hashtrue(real_lines)
	HashPredict,seq_dico=read_hashpredict(predicted_lines,label_seq_True)
	TP,FN,FP=counter(seq_dico,HashTrue,HashPredict,label_seq_True)
	Classiq_mesure(TP,FN,FP) 
	return 0			
#===================================================================================
#			    		Main
#===================================================================================
def main():
    usage = usage = "python sequence_based_fscore.py -p <clustering output> -t <true cluster file> \n"
    parser = OptionParser(usage)
    parser.add_option("-p", "--Predicted_Files_File", dest="Predicted_Files_File",
          help="read clusters from Predicted_File")
    parser.add_option("-t", "--True_clusters_file", dest="True_clusters_file",
          help="read data from True clusters file")
    
    (options, args) = parser.parse_args()
    if len(sys.argv) != 5:
        parser.error("incorrect number of arguments")
    
    Predicted_File = options.Predicted_Files_File
    True_file = options.True_clusters_file
    precision_recall(Predicted_File,True_file)

#===================================================================================
if __name__ == "__main__":
    main()
