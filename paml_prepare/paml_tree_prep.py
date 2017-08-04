
import inspect 
from Bio import Phylo
from Bio import Seq
from Bio import AlignIO
from Bio.AlignIO import  PhylipCodemlIO
from Bio.Alphabet import generic_dna
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
#from Bio import PhylipCodemlIO 



def message(case , variables=None):
    if case == "missing_tree_tip":
        seq , aln, tree = variables 
        message = """
        %s occurs in alignment file: %s
        but is not in %s:
        Check to make sure that your tree
        file is correct.
        \n""" %( seq , aln, tree )
    elif case == "help":
        message="""
        paml_tree_prep.py -f <fasta_files> -t <tree_file>
        
        outputs strict phyllip file and reduced tree file for use with PAML.
        """

    return message


def edit_for_paml_codon( seq_record):
    ret_seq =seq_record
    seq_len = len(seq_record.seq)
    seq_len = seq_len - (seq_len %3 )
    ret_seq.seq = seq_record.seq[0:seq_len]
    return ret_seq                     




def get_ghost_and_living(aln_file, aln_type):
    """
    ghost seqs are place holder seqs added in another program
    they consist of only "-" characters. Living seqs are seqs
    that contain aligned sites. This function returns a set of 
    seq names to include and a reduced alignment for analysis in
    PAML or other set.     
    """
    keep_seqs = []
    keep_seq_id = []
   
   
    for seq_record in list(AlignIO.parse(aln_file, aln_type))[-1]:
       gaps = str(seq_record.seq).count("-")
       seq_len = len( seq_record.seq)
       if seq_len != gaps:
           seq_record = edit_for_paml_codon(seq_record)
           keep_seqs.append(seq_record)
           keep_seq_id.append(seq_record.id)
    red_aln = MultipleSeqAlignment(keep_seqs)
    return [red_aln, keep_seq_id]


def collapse_redundant_nodes(tree):
    """
    function takes Tree object, tree, and counts number
    of tips per clade, if the number of 
    tips in the current clade == the 
    number of tips in the previous clade
    the clade is collapsed because it is
    redunant. Rendundant clades can arise
    because of user error or as in the case
    of reduce_tree(), all the tips are removed
    from a clade during a trimming process. 

    """
    counter = 0
    prev_count = None
    remove_list = []
    clade_dict = {}
    for clade in tree.get_nonterminals():
        print (counter, type(counter))
        count = clade.count_terminals()
        clade_dict.setdefault( count, []).append(clade)

    for count_key in clade_dict.keys():
        if len(clade_dict[count_key]):
            for subclade in clade_dict[count_key]:
                for clade in clade_dict[count_key]:
                    if subclade in clade and subclade not in remove_list:
                        remove_list.append(subclade)

    for sub_clade in remove_list :
        tree.collapse(sub_clade)
    
    return tree
                        
            




def reduce_tree(tree,keep_seq_id):
    """
    reduce tree takes a list of all alignment names
    that occur within a tree. It elminates any names
    not on that list, creating a tree that reflects the
    membership of an alignment. 
    """
    added_list=[]
    for terminal_name  in tree.get_terminals():
        if str(terminal_name) not in keep_seq_id:
            tree.prune(terminal_name)
        else :
            added_list.append(str(terminal_name))
    print(tree)
    for seq in keep_seq_id:
        assert( seq in added_list ), message("missing_tree_tip" ,[ seq, aln_file, tree_file ])

    collapse_redundant_nodes(tree)
    Phylo.draw_ascii(tree)
    return tree

 
def produce_phy_tree(aln_file, aln_type, out, exon, tree_type, tree):
    """
    run all function to create an inframe exon specific paml_file
    or create an intron (non-coding) phylip file.
    """
    
    red_aln, keep_seq_id  = get_ghost_and_living(aln_file,aln_type)
    if exon is True:
        AlignIO.write( red_aln,  "%s.phy" %(out), "phylip-codeml-g" )
    else:
       AlignIO.write( red_aln,  "%s.phy" %(out), "phylip-codeml" ) 

    Phylo.write(reduce_tree(tree,keep_seq_id), "%s.tre" %(out), tree_type) 
    Phylo.draw_ascii(tree)


if __name__ == "__main__" :
    import argparse
    import copy
    parser = argparse.ArgumentParser(description=message("help"))
    parser.add_argument('-f', '--fasta',
               nargs='+',
               help="space delimited list of aligned fasta files", 
               type=str,
               dest="fasta_files")
    parser.add_argument('-t', '--tree_file',  
                help="Newick tree file that contains all possible taxons in the alignment in essence this is a master a file that we reduce to contain only the taxa that occur for each tree.", 
                type=str,
                dest='tree_file')
    parser.add_argument('-o', '--out', 
                help='out prefix for fasta file[s] in list we create a new PAML compatible phyllip file, and newick tree file..',
                type=str,
                dest='out',
                default="paml_prep_")
    parser.add_argument('--type',
                        help=' type of input file to open. Default is fasta, and availble types are limited to those available in biopython 1.7 in python 3.',
                        type=str,
                        default='fasta',
                        dest="aln_type")
                          
    parser.add_argument('--exon',
                        help='use flag if sequences are in frame and you wish to have a codon position aware paml output file.\nPlease note that this depends on a modified version of biopython module provided in virtual environment.',
                        action='store_true')
    parser.add_argument('--tree_type',
                        help='It is possible to use other tree file formats amenable to biopython 1.7.',
                        default="newick" )


    argv = parser.parse_args()

    


    tree_file= argv.tree_file
    tree_type = argv.tree_type
    out = argv.out
    aln_type = argv.aln_type
    steady_tree = Phylo.read( argv.tree_file, argv.tree_type)
    
    if argv.exon :
        exon = True
    else:
        exon = False


    for fasta in argv.fasta_files :
        tree = copy.deepcopy(steady_tree)
        out = "%s%s" %(argv.out , ".".join(fasta.split("/")[-1].split(".")[0:-1]))
        produce_phy_tree(fasta, aln_type, out, exon, tree_type, tree)

    

    Phylo.draw_ascii(tree)
       


