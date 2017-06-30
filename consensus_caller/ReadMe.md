**pysam_consensus.py** is a simple consensus caller it outputs the most abundant characters, base, insertion, or deletion for each mapped nucleotide position. If the most abundant base falls below preset cutoff default 20, base is reported as character N. Sequences are padded to the first mapped read position in the contig with character P. There is no internal padding for bases with zero coverage. These are not reported and can look like deletions which are not reported, since they are not present in the sequence. As with all consensus calling from maps, the quality of consensus sequence is reliant on the quality of the map from which it is derived. It is advisable to remove duplicates, and realign indels using tools such as picard and GATK prior to calling a consensus. 


**requires pysam**
** written in python3.4 venv by NDHall**




```
    usage: pysam_consensus.py -b <bam_file>[required] -o <ouput>[required] -d <depth>[default=20] -h <--help> this message

     pysam_consensus.py calls a consensus by iterating over all the pileup columns in an indexed bam file.
     default minimum depth can be set from the command line via -d or --depth flag. Each sequence is padded
     with character P to the first position at which mapped reads occur in the pileup. This keeps sequences
     in frame and can easily be turned into spaces or deleted via sed '/>/! s/P/-/g' or sed '/>/! s/P//g'. 

     Flags

    -h --help     this message
    
    -b --bam_file indexed bam file.
   
    -o --output   output file as you intended it to be named. 

    -d --depth    minimum acceptable depth of reads from whic
                  to call a consensus. For position 1 of fasta
                  to be called A there must be 20 A's recovered
                  from that column. If A is the the number with 
                  the maximum characters and it occurs 19 times
                  position 1 will be set as N. 


```