
import sys
import tqdm
import math
import time
import resource
import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
from collections import Counter
from optparse import OptionParser

####################################################################
"""read the fasta file of input sequences"""	
def readFastaMul(nomFi):
	f=open(nomFi,"r")
	lines=f.readlines()
	f.close()
	seq=""
	nom=""
	lesSeq={}
	for l in lines:
		if l[0] == '>':
			if seq != "":
				lesSeq[nom] = seq
			nom=l.split("\n")[0][1:]
			seq=""
		else:
			seq=seq+l[:-1]
	if seq != "":
		lesSeq[nom.rstrip()] = seq.rstrip()
	return lesSeq

#####################################################################
# read the clustering result
def readClusteringResults(nomFi):
	# read the clustering result
	f=open(nomFi,"r")
	lines=f.readlines()
	f.close()
	nom=""
	Clustering_lables={}
	for l in range(len(lines)-1):
		cluster=lines[l].split("\t")[0].rstrip()
		Seq_nom=lines[l].split("\t")[1].split(" ")
		Seq_nom[-1]= Seq_nom[-1][:-1]
		Clustering_lables[cluster] = Seq_nom
	return Clustering_lables

#####################################################################
def gini(arr):
	## first sort
	sorted_arr = arr.copy()
	sorted_arr.sort()
	n = arr.size
	coef_ = 2. / n
	const_ = (n + 1.) / n
	weighted_sum = sum([(i+1)*yi for i, yi in enumerate(sorted_arr)])
	return coef_*weighted_sum/(sorted_arr.sum()) - const_

#####################################################################
def Plot(Clustering_lables,FastaFile):
	file = open(FastaFile.split(".")[0] +"_cluster_distribution.txt","w")
	fig = plt.gcf()
	fig.subplots_adjust(hspace=0.4)
	fig.set_size_inches((20,12))
	fig.suptitle(FastaFile.split(".")[0], fontsize=20)
	x=[]

	for key in Clustering_lables.keys():
		x.append(len(Clustering_lables[key]))
		if len(Clustering_lables[key]) == 0:
			print key
	X=np.array(sorted(x))


	ax = fig.add_subplot(221)
	N = len(x)
	y = [0]*len(x)
	s = [(n/min(X)) for n in X]
	colors = np.random.rand(N)
	plt.scatter(X,y,s=s,c=colors,alpha=0.5)
	
	plt.axis('equal')
	plt.xlabel('sequence number in cluster')

	ax = fig.add_subplot(222)
	l=plt.plot(X,'k.')
	plt.setp(l, markersize=3)
	ax.set_yscale("log", nonposy='clip')
	plt.xlabel('Cluster')
	plt.ylabel('Log of sequence number in cluster')

	ax = fig.add_subplot(223, aspect='equal')
	X_lorenz = X.cumsum() /float( X.sum())
	X_lorenz = np.insert(X_lorenz, 0, 0)
	X_lorenz[0], X_lorenz[-1]
	ax.scatter(np.arange(0,1,1.0/X_lorenz.size),X_lorenz, marker='.', color='darkgreen', s=100)
	ax.plot([0,1], [0,1], color='k')
	ax.text(0.0, 0.9, "Gini coefficient = %.2f" % gini(X),fontsize=10)
	plt.xlabel('Cumulative share of clusters')
	plt.ylabel('Cumulative share of sequences')

	ax = fig.add_subplot(224)
	rev=np.fliplr([X])[0]
	clone_list = list(rev)
	seq_number = sum(clone_list)
	j = 1
	for i in clone_list :
		line_to_write = "Cluster " + str(j) + "\t"+ str(i/float(seq_number)) + "\n"
		file.write(line_to_write)
		j += 1
	file.close()

	if len(X) > 100 :
		axe=list(range(1,101))
		s = [(float(n)/(sum(X))) for n in rev[0:100]]
		plt.bar(axe ,s)
		plt.axhline(y=0.05 ,linewidth=1, color='r')
		plt.xlim(0,100)
	else :
		axe=list(range(1,len(X)+1))
		s = [(float(n)/(sum(X))) for n in rev[0:len(X)]]
		plt.bar(axe ,s)
		plt.axhline(y=0.05 ,linewidth=1, color='r')
		plt.xlim(0,len(X))

	plt.xlabel('Cluster')
	plt.ylabel('relative Cluster size')
	fig.savefig(FastaFile.split(".")[0]+'.png', dpi=500)

#	plt.show()
#####################################################################
def main():
	usage = "usage: %prog -f FastaFile -c ClusteringFile"
	parser = OptionParser(usage)
	parser.add_option("-f", "--FastaFile", dest="FastaFile",
	      help="read data from FILENAME")
	parser.add_option("-c", "--ClusteringFile",dest="ClusteringFile",
	      help="read data from ClusteringFile")
	(options, args) = parser.parse_args()
	if len(sys.argv) != 5:
		parser.error("incorrect number of arguments")

	FastaFile = options.FastaFile
  	ClusteringFile = options.ClusteringFile
	Dicofasta=readFastaMul(FastaFile)
	Dicoresult=readClusteringResults(ClusteringFile)

	Plot(Dicoresult,FastaFile)
	print "Done!"

#####################################################################
if __name__ == "__main__":
	main()


