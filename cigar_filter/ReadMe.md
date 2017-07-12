pysam_cigar_filter.py was written in python 3.3 using pysam.
``` 
    pysam_cigar_filter.py  -s or --sam_file <sam or bam file> 
                           -l --length <Minimum acceptable match length int default=20> 
                           -o --out_prefix <out_prefix> 
                           -h --help [returns help option]
```

The program keeps,soft clipped reads that over lay contig ends and reads with soft clipping at one end and fewer
than 5 bases soft clipped. It filters out the rest of the soft clipped reads.

