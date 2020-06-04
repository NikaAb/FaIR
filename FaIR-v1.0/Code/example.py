import time
import sys
from clusterIO import readFasta, write_cluster, print_unmerged
from ClusterCDR3Only import MixClust_CDR3Only



def example_run():
    #print args
    D,M,O,I,MI,T,SS=("example.fasta","CDR3","example_output",95.0,95.0,0,300)
    data=readFasta(D)
    print ("Total Sequences: ", len(data),'\n')
    if M=="CDR3":
        t=time.localtime()
        tm=(t.tm_hour, t.tm_min, t.tm_sec)
        
        t1 = time.time()
        #print(t1)
        clusters = MixClust_CDR3Only(data,th=I, tolerance=T, mth=MI, split_size=SS)

        #t=time.localtime()
        

        #print("Starting time: ",t.tm_hour, t.tm_min, t.tm_sec)
        t3=(time.time() - t1)/60.0
        
        print ("=========== Clustering Details==========\n")
        write_cluster(clusters, O+"_id_"+str(I)+".txt")
        print ("Total #Sequences: ", len(data),'\n')
        print ("Total #Clusters: ",len(clusters), '\n')
        
        print ("Clusters are written in <",  O+"_id_"+str(I)+".txt", "> file",'\n')
        t=time.localtime()
        print("Clustering started at ", tm[0],':',tm[1],':', tm[2], '\n')
        print("Clustering ends at ", t.tm_hour,':', t.tm_min,':',t.tm_sec, '\n')
        print ('Total Clustering Time:',round(t3,4) , "Mins\n")
        print("===========================================")
        print ("Thank you for using MixClust\n")


if __name__ == '__main__':
	example_run()