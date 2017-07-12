**fasta-stats.py** is modified from work by 
*Uli Koehler & Anton Smirnov copyright 2013*. I have modified it to produce a csv written to standard out that can then be filtered by **fasta_stats_parse.py**.

Usage: 
``
fasta_stats.py <fastas…>
``


**fasta_stat_parse.py** is designed to handle the output from print_fasta_stats.py a python script that tallies A,C,T,GN,P characters in fasta. It is designed to be used downstream from **pysam_consensus.py** and **fasta-stats.py**, which uses N to indicate bases in the original bam map that fell below minimum acceptable coverage and P to indicate bases with no read coverage in the orginal bam map. Note, P is not the same as deletion cigar. Deletion cigars are handled in the context of **pysam_consensus.py** and regions marked by P for padding are not overlapped by reads. 

Usage: 
```
  fasta_stat_parse.py -c <csv>[required] --max_N_rel_freq  <default:4>   \ 
                                           --max_N <default:9>              \ 
                                           --max_P_rel_freq <default:1> 
        
       -c --csv                       csv file prodced by use of print_fasta_stat.py 
                                     on pysam_consensus.py output
       
       --max_N_rel_freq = default:4  where N is number of bases that did not meet 
                                     coverage requirement
       
       --max_N          = default:9  hard limit for N free of propotional expansion 
                                     for very larger sequences
       
       --max_P_rel_freq = default:1  where P is number of padded bases that have no 
                                     coverage. 
                                     
```                                     
Long spans of N or P can make alignment of sequences problematic in clustal and muscle run through megacc. By limiting number we can maximize the probability that alignments are correct
