##get_insert.py 

get_insert.py is a set of functions that uses pysam to emulate mpileup.
Since pysam does not support insertion reporting in its pileup emulator
an extra function is needed to get insertion sites if they exist.

Briefly, as pysam iterates over a pileup by column, attribute indel is 
used to check if an insertion is next. If an insertion is next, the cigar
tuples for the read are called, parsed and used to extract the insertion 
sequence from the raw read. A single 5'base +insertion region are reported.
reported region is added to list which is counted using a Counter function.

Deletions are detectable because the read in pileupcolumn has attribute .is_del.
These are noted as D for the purposes of counting and removed prior to writing fasta. Each seperate contig in the map is padded to the first starting position with character P. Uncovered based in the middle or end of pileup are not padded. This could be added later. Bases that are present, but that fall below minium cutoff value are reported as character N.  

##10K.reps

When breaking ties, among equally supported bases for the consensus, get_insert.py relies on the unordered and pseudo-random nature of python dictionaries. 10K.reps is a fasta from 1 map created 10K times. With multiple equally supported bases in the middle. sort| uniq -c of sequence reveals ~ equal distribution of all possible bases.
```
   
   2388 TAAAAAAAAAAAAAAAAAAAATTTTTTTTTTTTTTTTTTTTA
   2417 TAAAAAAAAAAAAAAAAAAAACTTTTTTTTTTTTTTTTTTTTA
   2572 TAAAAAAAAAAAAAAAAAAAACATTTTTTTTTTTTTTTTTTTTA
   2623 TAAAAAAAAAAAAAAAAAAAACGTTTTTTTTTTTTTTTTTTTTA

```



## bam file
bam file was used to make 10K.reps

## test.fasta
is test output fasta
