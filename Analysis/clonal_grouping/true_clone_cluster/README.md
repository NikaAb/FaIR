# Silhouette
**Command-line manual**







# Input 

The input files of Silhouette.py are: 

Clustering_File which contains the output of clustring tool. This file should have the folowing format :

``` diff
Cluster_number	seq1 seq2 seq3,...

Example :

1	S483 S483 S483
```

FastaFile which contains nucleotide sequences. This file should have the folowing format :

``` diff
> seq Id
seq

Example :

>S953
CAGGTCCAGCTTGCGCAGTCTGGGGCTGAGGTGAAAAAGCCTGGGGCCTCGGTGAAGGTTTCCTGCAAGGCTTCTGGATACACCTCCAC
```


# Output

## Main output 

The output of Silhouette.py is the value of silhouette:

``` diff
Silhouette : 0.68  
```

## Usage

``` bash
$ python Silhouette.py -f Patient_one.fa  -c Patient_one_Clone_formatted.txt 

   ```

