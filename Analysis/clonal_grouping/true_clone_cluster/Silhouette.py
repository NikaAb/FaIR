"""
Evaluating the quality of any clustering tool even on sequences that we dont know the real cluster labeles.
"""
#Have to install optparse and tqdm

import sys
import tqdm
import math
import time
import resource
import Levenshtein
import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt

from collections import Counter
from optparse import OptionParser
from sklearn.cluster import KMeans
from sklearn.cluster import AffinityPropagation


#####################################################################
def readFastaMul(nomFi):
	"""read the fasta file of input sequences"""	
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
			nom=l[1:-1]
			seq=""
		else:
			seq=seq+l[:-1]
	if seq != "":

		lesSeq[nom.rstrip()] = seq.rstrip()

	return lesSeq

#####################################################################

# read the clustering result

def readClusteringResults(nomFi):
	f=open(nomFi,"r")
	lines=f.readlines()
	f.close()
	Clustering_lables={}
	for l in range(len(lines)):
		Seq_nom=lines[l].split("\t")[1].rstrip().split(" ")
		cluster=lines[l].split("\t")[0].rstrip()

		if cluster in Clustering_lables.keys():
			Clustering_lables[cluster].append(Seq_nom)
		else:
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
	
	#fig = plt.figure()
	fig = plt.gcf()
	fig.subplots_adjust(hspace=0.4)
	fig.set_size_inches((20,12))
	fig.suptitle(r' Simulated repertoir '+FastaFile.split(".")[0], fontsize=20)
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
	plt.xlabel('Number of sequence in cluster')

	
	
	ax = fig.add_subplot(222)
	l=plt.plot(X,'k.')
	plt.setp(l, markersize=3)
	ax.set_yscale("log", nonposy='clip')
	plt.xlabel('Number of clusters')
	plt.ylabel('Log of number of sequence in cluster')



	ax = fig.add_subplot(223, aspect='equal')
	X_lorenz = X.cumsum() /float( X.sum())
	X_lorenz = np.insert(X_lorenz, 0, 0)
	X_lorenz[0], X_lorenz[-1]
	ax.scatter(np.arange(0,1,1.0/X_lorenz.size),X_lorenz, marker='.', color='darkgreen', s=100)
	ax.plot([0,1], [0,1], color='k')
	ax.text(0.0, 0.9, "Gini coefficient = %.2f" % gini(X),fontsize=10)
	plt.xlabel('Cumulative share of clusters')
	plt.ylabel('Cumulative share of sequences')
	
	fig.savefig(FastaFile+'.png', dpi=500)
	

	#plt.show()
	
#####################################################################
def CalculateMedoid(Dicofasta,Dicoresult):
	centroid={}
	print "Calculating Medoid sequence of each cluster ... \n"
	for key in tqdm.tqdm(Dicoresult.keys()) :
		listloc=[]
		for seq in Dicoresult[key]:
			if seq.rstrip() in Dicofasta.keys():
				listloc.append(Dicofasta[seq.rstrip()])
			else :
				print "Caution the",seq.rstrip(),"is not in the fasta file"
				
		centroid[key]=Levenshtein.median(listloc)

	return centroid
#####################################################################
def CalculateMedianDist(Dicocentroid):
	Cluster_list=[]
	Centroid_list=[]
	dicoNeighbour={}	
	for cluster in Dicocentroid.keys():
		Cluster_list.append(cluster)
		
		Centroid_list.append(Dicocentroid[cluster])
	print "Creating neighnerhood dictionary ... \n"
	list__for_final_elem=[]
	for i in tqdm.tqdm(range(len(Centroid_list)-1)) :
	#for i in range(len(Centroid_list)-1):
		#print i
		listloc=[]
		for j in range(i+1,len(Centroid_list)):
			
			listloc.append(Levenshtein.distance(Centroid_list[i],Centroid_list[j]))
			
			if j == len(Centroid_list)-1:
				list__for_final_elem.append(Levenshtein.distance(Centroid_list[i],Centroid_list[j]))

		argmin=listloc.index(min(listloc))
		if Cluster_list[i] not in dicoNeighbour.keys():
			dicoNeighbour[Cluster_list[i]]=Cluster_list[argmin+i+1]
		if Cluster_list[argmin+i] not in dicoNeighbour.keys():
			dicoNeighbour[Cluster_list[argmin+i+1]]=Cluster_list[i]
	dicoNeighbour[Cluster_list[-1]]=Cluster_list[list__for_final_elem.index(min(list__for_final_elem))]

	return dicoNeighbour

#####################################################################

# SILHOUETTE

def silhouette(Dicofasta,Dicocentroid,Dicoresult,DicoNeighbour):
	
	summe=0
	print "Calculating silhouette ... \n"
	for cluster in tqdm.tqdm(Dicoresult.keys()) :
	#for cluster in Dicoresult.keys():
		for seq in Dicoresult[cluster]:
			ai=Levenshtein.distance(Dicofasta[seq],Dicocentroid[cluster])
			bi=Levenshtein.distance(Dicofasta[seq],Dicocentroid[DicoNeighbour[cluster]])
			summe+=calculeSil(ai,bi)
 	return summe/float(len(Dicofasta))

#####################################################################

def calculeSil(ai,bi):
    si=0
    if ai!=bi:
        if ai<bi:
            si=1-(ai/float(bi))
        else:
            si=(float(bi)/ai)-1
    return si


####################################################################
def main():
	usage = "usage: Silhouette.py -f FastaFile -c ClusteringFile"
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
	time_start = time.clock()
	Dicofasta=readFastaMul(FastaFile)
	#print Dicofasta

	Dicoresult=readClusteringResults(ClusteringFile)


	#Plot(Dicoresult,FastaFile)
	
	Dicocentroid=CalculateMedoid(Dicofasta,Dicoresult)

	DicoNeighbour= CalculateMedianDist(Dicocentroid)

	sil=silhouette(Dicofasta,Dicocentroid,Dicoresult,DicoNeighbour)
	print "Silhouette :", sil

	time_elapsed = (time.clock() - time_start)

	print "Calculation time : ",time_elapsed

#####################################################################
if __name__ == "__main__":
	main()
  


