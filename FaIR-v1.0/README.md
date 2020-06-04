# FaIR 2019.01

**Command-line manual**

## About
FaIR is an ultra-fast and sensitive algorithm for clustering BCR-IG sequences. Based on heuristic yet high-performance methods, FaIR can cluster BCR-IG sequences before VDJ assignment in a reasonable time.

# Requirements and installation
The Levenshtein Python C extension package, you can install it by using pip :

``` bash
 pip install python-Levenshtein
```


# Input and parameters

The main input file of FaIR is a *set of reads*, given as a `.fasta`.


## Main algorithm parameters
``` diff
  -dataset                 Name of input fasta file.
  -output                  Name of output file.
  -method                  Clustering method (default: CDR3, based on CDR3 sequence).
  -id                      Similarity threshold
  -splitsize               Size of a smaller block of sequences. It helps to speed up the clustering by addressing a chunk of sequence to each available core (default: 300)
```


# Output

## Main output files

The main output of FaIR is a three column text file, seperated by space, in the following format:
``` diff
  Cluster_number  sequence_id sequence 
  1 S11_1 GTTTTTCTGTTCACAGGGGTCCTGTCCCCGGTGCAGCTACAGCAGTGGGGCGCAGGACTGTTGAAGCCTTCGGAGACCCTGTCCCTCACCTTCGCTGTCTATGGTGGGTCCTTCAGTGGTTACTACTGGAGCTGGTTCCGCCAGCCCCCAGGGAAGGGGCTGGAGTGGATTGGGGAAATCAATCATAGTGGAAGCACCAACTACAACCCGTCCCTCAAGAGTCGAGTCACCATATCAGTAGACACGTCCAAGACC
```

## Example

``` diff
& python example.py
```


## Usage

``` bash
$ python FaIRt.py -dataset simulated_seq.fa -output simulated_seq_out  -method CDR3 -id 90.0 -splitsize 300
   # clustering based on junction region (which contains CDR3)
   # Similarity threshold for creating initial clusters is 90%.
   # The result of this clustering is in the simulated_seq_out_id_90.txt file ('clusters')
   # sequences are divided into groups of 300, and each group is given to one core.
   ```
