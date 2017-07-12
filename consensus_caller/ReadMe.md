**pysam_consensus.py** is a simple consensus caller. It outputs the most abundant characters,( base, insertion, or deletion ) for each mapped nucleotide position. If the most abundant base falls below the preset cutoff, default 20, the base is reported as character N. Sequences are padded from the begining of the contig to the first mapped read position in the contig with character P and there is internal padding for regions that are not covered by any reads. 
 


**requires pysam
written in python3.4 venv by NDHall**




```
    usage: pysam_consensus.py -b <bam_file>[required] -o <ouput>[required] -d <depth>[default=20] -h <--help> this message

     pysam_consensus.py calls a consensus by iterating over all the pileup columns in an indexed bam file.
     default minimum depth can be set from the command line via -d or --depth flag. Each sequence is padded
     with character P to the first position at which mapped reads occur in the pileup. This keeps sequences
     in frame and can easily be turned into spaces or deleted via sed '/>/! s/P/-/g' or sed '/>/! s/P//g'. 

     Flags

    -h --help     this message
    
    -b --bam_file indexed bam file
   
    -o --output   full output name 

    -d --depth    minimum acceptable depth of reads from which
                  to call a consensus. For position 1 of fasta
                  to be called A there must be 20 A's recovered
                  from that column. If A is the the number with 
                  the maximum characters and it occurs 19 times
                  position 1 will be set as N. 


```
