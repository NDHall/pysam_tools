
##fasta_ghost.py



fasta_ghost.py is to prepare a series of fasta files for 
multiple alignments. It takes a list file, one or more 
fasta files and an outprefix. It reads each fasta and
if a fasta is missing an accession a new blank accession
of the appropriate length is added. Accessions are added in
the order they are listed in the SEQ_LIST_FILE. 

usage:
```
fasta_ghost.py [-h] [-f FASTA_FILES [FASTA_FILES ...]]
                      [-l SEQ_LIST_FILE] [-o OUT]

```


