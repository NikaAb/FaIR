# Junction region extraction 

In order to extract the CDR3 region from each sequence, we adapt the algorithm 
proposed in vidjil software. 

Vidjil builds a hash-table with a set of k-words extracted from IMGT-gapped V and
J germline sequences. The choice of k depends on the germline to be indexed, and
it is chosen to minimise the number of shared words between V and J genes.The hash
-table is then used to index  IGH input sequence by looking for the boundary of 
V and J genes. 

Trimming the IGH sequence closes to the VJ boundaries, theoretically, gives us the
CDR3 region, which contains the N1-D-N2 sequence. However, D genes usually have a
very short length, and the extracted segment might not contain enough information,
that could  compromise the clustering performance. To overcome this problem, we 
first tried to extend the extracted segment found by vidjil. But, it requires to
automatically define the length of a such extension for both directions. Some tests
have shown that infers these parameters could burden the run-time of our approach.
Then, we have modified the germline indexing by trimming the V and J germline genes,
differently.

For V germline genes we keep only the segment before the Cysteine at position 104
in FR3-IMGT region, and for J genes we keep only the segment after the Phenylalanine
or Tryptophane at position 118 in FR4-IMGT region. By using these shorter germline 
sequences, we can assure that the extracted segment will cover at least the CDR3 
sequence.

You can find the modified germline in the vidjil_master_modified folder.
***

## 1- Running vidjil algo with the trimmed germline
```
$ time vidjil_master_modified/vidjil-algo -c windows -U  -g germline/homo-sapiens.g monoclonal.fasta
real    0m4.791s
user    0m3.881s
sys    0m0.626s
```
## 2- Extracting CDR3 sequences from the vidjil algo's output 
```
$ time python junction_extraction.py -v monoclonal.segmented.vdj.fa -f monoclonal.fasta -o monoclonal.cdr3.fa
real    0m1.229s
user    0m0.441s
sys    0m0.407s
```
