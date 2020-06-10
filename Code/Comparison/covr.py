import sys
import tqdm
import math
import time
import resource
import Levenshtein
import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt


import collections 
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
# read the clustering result
def readClusteringResults(nomFi):
	lines = read_file (nomFi)
	Clustering_lables, seq_clust_corresp = {}, {} 
	for l in range(len(lines)):
		Seq_nom = lines[l].split("\t")[1].rstrip().split(" ")
		cluster = lines[l].split("\t")[0].rstrip()
		for seq in Seq_nom :
			seq_clust_corresp[seq] = cluster
		
		if cluster in Clustering_lables.keys():
			Clustering_lables[cluster].append(Seq_nom)

		else:
			Clustering_lables[cluster] = Seq_nom
	return Clustering_lables, seq_clust_corresp
#####################################################################
def remove_empty_keys(d):
    for k in d.keys():
    	if d[k] == [] :
            del d[k]
#####################################################################
# label all the sequences in the second distribution with join, split, identical
def label_seq(dico_dstb_one, dico_dstb_two, seq_clust_corresp_one, seq_clust_corresp_two):
	dico_dstb_two_seq = {}
	Mix_cluster_detail = {}
	not_found = []
	for key in dico_dstb_two.keys():
		#print key,"cluster in file 2" 
		Mix_cluster_detail[key] = []
		list_loc_join = []
		list_loc_split = []
		for seq in dico_dstb_two[key] :
			if seq in seq_clust_corresp_one.keys() :
				corresponding_cluster_dstb_one = seq_clust_corresp_one[seq]
			else :
				not_found.append(seq)
			# find the identical
			if collections.Counter(dico_dstb_two[key] ) == collections.Counter(dico_dstb_one[corresponding_cluster_dstb_one]) :
				#print dico_dstb_two[key],dico_dstb_one[corresponding_cluster_dstb_one]
				dico_dstb_two_seq[seq] = 'identical'
			elif len(collections.Counter(dico_dstb_two[key] ) ) > len(collections.Counter(dico_dstb_one[corresponding_cluster_dstb_one]) ):
				#print dico_dstb_two[key],dico_dstb_one[corresponding_cluster_dstb_one]
				for s in  dico_dstb_one[corresponding_cluster_dstb_one] :
					list_loc_join.append(s)
			else :
				#print dico_dstb_two[key],dico_dstb_one[corresponding_cluster_dstb_one]
				list_loc_split.append(seq)
		if list_loc_join != []:
			if collections.Counter(dico_dstb_two[key] ) == collections.Counter(list_loc_join) :
				#print dico_dstb_two[key],list_loc_join
				for j in dico_dstb_two[key] :
					dico_dstb_two_seq[j] = 'join'
			else :
				# mix clusters (identical)
				list_mix = []
				for mix in dico_dstb_two[key] :
					#print mix,"mix"
					if mix in seq_clust_corresp_one.keys():
						list_mix.append(seq_clust_corresp_one[mix])
					else: 
						not_found.append(seq)
				#print "list_mix",collections.Counter(list_mix)
				for key_c in collections.Counter(list_mix).keys():
					#print key,seq_clust_corresp_one
					if len(dico_dstb_one[key_c]) == collections.Counter(list_mix)[key_c]:
						Mix_cluster_detail[key].append(key_c)
						for iden in dico_dstb_one[key_c] :
							dico_dstb_two_seq[iden] = 'identical'
					else :
						for spl in dico_dstb_one[key_c] :
							dico_dstb_two_seq[spl] = 'split'

		if list_loc_split != []:
			for l in dico_dstb_two[key] :
					dico_dstb_two_seq[l] = 'split'
	print "not founds : ", not_found
	remove_empty_keys(Mix_cluster_detail)
	return dico_dstb_two_seq ,Mix_cluster_detail
####################################################################
def label_cluster(dico_dstb,dico_dstb_seq,Mix_cluster_detail,first_distribution,second_distribution):
	cluster_label_list = []
	cluster_label_count = {"Mix" : 0 , "split" : 0, "identical" : 0, "join" : 0}
	#print Mix_cluster_detail
	for key in dico_dstb.keys() :
		list_loc = []
		for seq in dico_dstb[key] :
			if seq in dico_dstb_seq.keys():
				list_loc.append(dico_dstb_seq[seq])
			else: 
				print seq
		if len(list(set(list_loc))) == 1:
			cluster_label_list.append(list_loc[0])
			cluster_label_count[list_loc[0]] += len(dico_dstb[key])
		else :
			#print list_loc
			cluster_label_list.append("Mix")
			cluster_label_count["Mix"] += len(dico_dstb[key])
	#print cluster_label_count, 
	label_dico = dict(collections.Counter(cluster_label_list))
	plt.bar(range(len(label_dico)), label_dico.values(), align='center')  # python 2.x
	plt.xticks(range(len(label_dico)), label_dico.keys())  # in python 2.x
	plt.ylabel('# cluster')
	#plt.ylim((0, 170)) 
	#cluster_one = first_distribution.split("_")[0]  + "_" + first_distribution.split("_")[3].split(".")[0] # simulation 
	#cluster_two = second_distribution.split("_")[0] + "_" + second_distribution.split("_")[3]  # simulation
	cluster_one = "I1 IMGT"#first_distribution.split(".")[0]
	cluster_two = "FaIR"#second_distribution.split(".")[0]
	name = "Distribution comparison : "+ str(cluster_one)+ " and "  +str(cluster_two)
	plt.title(name)
	annotation = []
	for k in label_dico.keys() :
		#print k
		annotation.append("{:.2f}".format(float(cluster_label_count[k])/sum(cluster_label_count.values())))
	for x,y,z in zip(range(len(label_dico)),label_dico.values(),annotation):
		#print x,y
		plt.annotate(z, # this is the text
	                 (x,y), # this is the point to label
	                 textcoords="offset points", # how to position the text
	                 xytext=(0,10), # distance from text to points (x,y)
	                 ha='center') # horizontal alignment can be left, right or center

	plt.show()
	plt.save()
	"""
	fig, ax = plt.subplots(figsize=(6, 3), subplot_kw=dict(aspect="equal"))

	recipe = ["225 g flour",
	          "90 g sugar",
	          "1 egg",
	          "60 g butter",
	          "100 ml milk",
	          "1/2 package of yeast"]

	data = [225, 90, 50, 60, 100, 5]

	wedges, texts = ax.pie(data, wedgeprops=dict(width=0.5), startangle=-40)

	bbox_props = dict(boxstyle="square,pad=0.3", fc="w", ec="k", lw=0.72)
	kw = dict(arrowprops=dict(arrowstyle="-"),
	          bbox=bbox_props, zorder=0, va="center")

	for i, p in enumerate(wedges):
	    ang = (p.theta2 - p.theta1)/2. + p.theta1
	    y = np.sin(np.deg2rad(ang))
	    x = np.cos(np.deg2rad(ang))
	    horizontalalignment = {-1: "right", 1: "left"}[int(np.sign(x))]
	    connectionstyle = "angle,angleA=0,angleB={}".format(ang)
	    kw["arrowprops"].update({"connectionstyle": connectionstyle})
	    ax.annotate(recipe[i], xy=(x, y), xytext=(1.35*np.sign(x), 1.4*y),
	                horizontalalignment=horizontalalignment, **kw)

	ax.set_title("Matplotlib bakery: A donut")

	plt.show()
	"""
####################################################################
def main():
	usage = "usage: covr.py -f first_distribution -g second_distribution"
	parser = OptionParser(usage)
	parser.add_option("-f", "--first_distribution", dest="first_distribution",
	      help="read data from first_distribution")
	parser.add_option("-g", "--second_distribution",dest="second_distribution",
	      help="read data from second_distribution")
	(options, args) = parser.parse_args()
	if len(sys.argv) != 5:
		parser.error("incorrect number of arguments")
	
	first_distribution = options.first_distribution
	second_distribution = options.second_distribution
	one_Clust_lables, seq_clust_corresp_one = readClusteringResults(first_distribution)
	two_Clust_lables, seq_clust_corresp_two = readClusteringResults(second_distribution)
	#print one_Clust_lables, seq_clust_corresp_one
	dico_dstb_seq ,Mix_cluster_detail = label_seq(one_Clust_lables, two_Clust_lables, seq_clust_corresp_one, seq_clust_corresp_two)
	label_cluster(two_Clust_lables,dico_dstb_seq,Mix_cluster_detail,first_distribution,second_distribution)
#####################################################################
if __name__ == "__main__":
	main()
  
