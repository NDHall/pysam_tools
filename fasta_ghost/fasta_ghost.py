#!~/Documents/pysam_consensus/venv/bin/python3


import argparse
from Bio import SeqIO






def message(case, variables=None):
    if case == "duplicate_accessions":
        to_add, seq_list_file
        message="""
     $s is a duplicate accession in
     file : %s 
     please remove duplicate accession.\n\n""" %( to_add, seq_list_file)

    elif case == "unequal_length" :
        rec, fasta = variables 
        message = """
     Check to make sure the file is aligned  
        %s does not have the same length as
        other sequences in the alignment 
        file: %s
       .
        \n""" %( rec , fasta)
    elif case == "help":
        message = """
        fasta_ghost.py is to prepare a series of fasta files for 
        multiple alignments. It takes a list file, one or more 
        fasta files and an outprefix. It reads each fasta and
        if a fasta is missing an accession a new blank accession
        of the appropriate length is added. 
        \n"""
    
    return message 





def get_seq_list(seq_list_file):
    f=open(seq_list_file, 'r')
    seq_split = f.read().split("\n")
    seq_list = []
    for accession in seq_split:
        to_add = accession.strip()
        if to_add != '':
            assert (to_add not in seq_list ), message("duplicate_accessions", [to_add, seq_list_file])
            seq_list.append(to_add)
    f.close()
    return (seq_list)


def get_seqs(fasta_file):

    """
    Function get_seqs() reads entire fasta into memory.
    If the sequences are very large and memory is limited
    this could be a problem.
    """

    seq_dict = {}
    prev_length = None 
    for rec in SeqIO.parse(fasta_file, 'fasta'):
        name = rec.id
        length = len(rec.seq)
        print(name, len(rec.seq))
        assert ( name not in seq_dict.keys() ), message("duplicate_accessions", [name, fasta_file])
        seq_dict[name] = rec
        if prev_length is not None:
            assert (length == prev_length ), message("uneqal_length",[rec.id, fasta_file]) 
    
        prev_length = length 
 
    return [length, seq_dict ]

    


def include_ghost_seqs( seq_list, seqs, length, out):
    f=open(out, "w")
    f.close()
    f=open(out, "a")
    
    for seq in seq_list :
        if seq not in seqs.keys():
            print(">%s\n%s" %( seq, "".join(["-"]*length)))
            f.write(">%s\n%s\n" %( seq, "".join(["-"]*length)))
        else:
            print( ">%s\n%s" %( seq, str( seqs[seq].seq)))
            f.write(">%s\n%s\n" %( seq, str( seqs[seq].seq)))

    f.close()
    



    


if __name__ =="__main__":
    parser = argparse.ArgumentParser(description=message("help"))
    parser.add_argument('-f', '--fasta',
               nargs='+',
               help="space delimited list of aligned fasta files", 
               type=str,
               dest="fasta_files")
    parser.add_argument('-l', '--list',  
                help="file containg lists of all accession names needed in the output fasta.", 
                type=str,
                dest='seq_list_file')
    parser.add_argument('-o', '--out', 
                help='out prefix for fasta file[s] in list.',
                type=str,
                dest='out',
                default="ghost_seqs")
    argv = parser.parse_args()

    seq_list = get_seq_list(argv.seq_list_file)
    print(  argv.fasta_files)
    for fasta in argv.fasta_files :
        
        out = "%s%s" %(argv.out , fasta.split("/")[-1])
        length, seq_dict = get_seqs(fasta)
        include_ghost_seqs( seq_list , seq_dict, length, out)

