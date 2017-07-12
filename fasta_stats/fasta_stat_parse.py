"""

For parsing output for characters by sequence per organism and determining
if organisms have enough high quaility sequences to continue in an analysis.
pandas_parse_fasta_stats.py output list of sequences with pass or fail indicated

"""

import numpy as np 
import pandas as pd
import sys
import getopt

def usage():
    message="""\n\n
    fasta_stat_parse.py is designed to handle the output from 
    print_fasta_stats.py a python script that tallies A,C,T,G
    N,P characters in fasta. It is designed to be used downstream
    from pysam_consensus.py, which uses N to indicate bases in 
    the original bam map that fell below minimum acceptable cov-
    erage and P to indicate bases with no read coverage in the
    orginal bam map. Note, P is not the same as deletion cigar.
    Deletion cigars are handled in the context of pysam_consensus.py
    and regions marked by P for padding are not overlapped by 
    reads. 

    Usage: 
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

    Long spans of N or P can make alignment of sequences problematic in clustal and 
    muscle run through megacc. By limiting number we can maximize the probability 
    that alignments are correct
    \n\n"""
    print ( message) 

def write_pass_fail(f_name,dataFrame):
    """
    requires dataFrame to contain 3 columns
    file_name = file name
    seq_name  = sequence name to be returned in files
    passed    = column indicating pass or fail of some
                required test.  
    """
    f_stem= ".".join(f_name.split(".")[:-1])+".list"
    pass_f_name = "pass_%s" %( f_stem )
    fail_f_name = "fail_%s" %( f_stem )
    dataFrame[(dataFrame.file_name == f_name)  & ( dataFrame.passed =='passed' )]['seq_name'].to_csv(pass_f_name,index=False   )
    dataFrame[(dataFrame.file_name == f_name)  & ( dataFrame.passed =='fail' )]['seq_name'].to_csv(fail_f_name,index=False   )



def call_pass_fail(fasta_stats, max_N_rel_freq=4, max_N = 9, max_P_rel_freq=1):
    """
    Column passed is determined based on 
        max_N_rel_freq = default:4  where N is number of bases that did not meet coverage requirement
        max_N          = default:9  hard limit for N free of propotional expansion for very larger sequences
        max_P_rel_freq = default:1  where P is number of padded bases that have no coverage. 
    Long spans of N or P can make alignment of sequences problematic in clustal and muscle run through
    megacc. By limiting number we can maximize the probability that our alignments are correct. 
  
    creating passed uses numpy.where(conditional , value_if_true, value_if_false )
    Each column tested in the pandas DataFrame has to be set off individually 
    for the test to work. numpy.where() supplants an overly complex application
    of apply.(lambda)

    Finally by calling it on a per sequence basis it is possible to now create best alignments for
    all the data we have available and then subsample, these down. 
    """

    fasta_stats["passed"] = np.where(((fasta_stats['N_rel_freq' ]<+ max_N_rel_freq) & (fasta_stats['N'] <= max_N ) & (fasta_stats['P_rel_freq'] <= max_P_rel_freq )), "passed","fail")
    ret_analysis = fasta_stats[["file_name","seq_name","passed"]]
    return ret_analysis


def split_by_name(ret_analysis):
    """ 
    Take list of file_name, seq_name, and passed.
    write out pass fail for each file in the parsed
    .csv file. using write_pass_fail():
    """

    for f_name in ret_analysis.file_name.unique():
        write_pass_fail(f_name, ret_analysis)



def main_argv_parse(argv):
    """
    argv_dict provides sensible defaults for 
        max_N
        max_N_rel_freq
        max_P_rel_freq
    which are the criteria for passing the 
    filtering step, these defaults are also 
    written into the function, but we can call 
    them here since this is where we interact 
    with command line. 
    
    """

    argv_dict={'csv':None , 'max_N_rel_freq':4, 'max_N':9, 'max_P_rel_freq':1};
    if len(argv[1:]) == 1 and ( argv[1] == '-h' or argv[1]== '--help'):
        usage()
        exit(0)
        
    try:
        opts, args = getopt.gnu_getopt(argv[1:],"h:c:",["help=","csv=","max_N_rel_freq=","max_P_rel_freq=","max_N="])

    except getopt.error:

        usage()
        sys.exit(2)
    for opt , arg in opts:
        if opt in ("-c","--csv"):
            argv_dict["csv"]=arg
        elif opt in ("--max_N_rel_freq"):
            argv_dict["max_N_rel_freq"]=int(arg)
        elif opt in ("--max_N"):
            argv_dict["max_N"]=int(arg)
        elif opt in ("--max_P_rel_freq"):
            argv_dict["max_P_rel_freq"]=int(arg)
        elif opt in ("-h","--help"):
            usage()
            exit(0)    
        else :
            print ("\n Option not recognized %s\n\n" %(opt))
            usage()
            
        
            sys.exit(1)
    if argv_dict['csv'] is not None:
        return argv_dict
    else:
        usage()
        exit(1)




#fasta_stats = pd.read_csv('/home/ndh0004/Dropbox/MitoGenes/Mega/GATK_V4_Ana/exon_fasta_stat.csv')

if __name__ == "__main__":
    argv_dict = main_argv_parse(sys.argv)
    fasta_stats = pd.read_csv(argv_dict["csv"])

    split_by_name(call_pass_fail(fasta_stats, argv_dict["max_N_rel_freq"],  argv_dict["max_N"],  argv_dict["max_P_rel_freq"] ))



