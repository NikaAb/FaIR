"""
Copyright (c) 2019 Bishnu Sarker (bishnukuet@gmail.com), Nika Abdollahi, Juliana Silva Bernardes

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""

"""
The command line instruction is:
python FaIR.py 
-dataset <Fasta file> 
-output <Prefix of output file> 
-method CDR3 
-id <similarity in % > 
-tolerance <residue Mismatch acceptable> 
-mid <Similarity Threshold to Merge clusters> 
-splitsize <Chunk Size>
Example:
python FaIR.py -dataset change_o_data.fa -output changeO_out -method CDR3 -id 90.0  -splitsize 300


"""
import time
import sys
from clusterIO import readFasta, write_cluster, print_unmerged
from ClusterCDR3Only import FaIR_CDR3Only



def get_option(arguments,option,pos=0):
    # Read the command line arguments

    l = len(arguments)
    i = pos
    while i<l:
        if arguments[i].upper() == option.upper():
            break
        else:
            i+=1
            continue
    return i

def parse_arguments(args):
    # Parse the command line argiments and retrieve the value
    k=get_option(args,"-output")
    output_prefix=args[k+1]
    k=get_option(args,"-dataset")
    dataset=args[k+1]
    k=get_option(args,"-method")
    method=args[k+1]
    k=get_option(args,"-id")
    id=float(args[k+1])
    k=get_option(args,'-tolerance')
    tolerance=int(args[k+1])
    k=get_option(args,'-mid')
    mid=float(args[k+1])
    k=get_option(args,'-splitsize')
    splitsize=int(args[k+1])


    return (dataset,method,output_prefix,id,mid,tolerance, splitsize)

def main():
    args=sys.argv[1:]
    #print args
    D,M,O,I,MI,T,SS=parse_arguments(args)
    data=readFasta(D)
    print ("Total Sequences: ", len(data),'\n')
    if M=="CDR3":
        t=time.localtime()
        tm=(t.tm_hour, t.tm_min, t.tm_sec)
        
        t1 = time.time()
        #print(t1)
        clusters = FaIR_CDR3Only(data,th=I, tolerance=T, mth=MI, split_size=SS)

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
        print ("Thank you for using FaIR\n")



if __name__ == '__main__':
    main()