

import get_insert # main functions live in get_insert
import getopt
import sys



def usage():
    message="""

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


    
    \n"""
    print ( message)

def main_argv_parse(argv):
    
    argv_dict={"output":'' ,"bam_file":'',"depth":20};
    if len(argv[1:]) == 1 and ( argv[1] == '-h' or argv[1]== '--help'):
        usage()
        exit(0)
        
    try:
        opts, args = getopt.gnu_getopt(argv[1:],"h:b:d:o:",["help=","output_prefix=","bam_file=","depth="])
    except getopt.error:
        print("\n\tFlag error see usage:\n")
        usage()
        sys.exit(2)
    for opt , arg in opts:
        if opt in ("-o","--output"):
            argv_dict["output"]=arg
        elif opt in ("-b","--bam_file"):
            argv_dict["bam_file"]=arg
        elif opt in ("-d","--depth"):
            argv_dict["depth"]=int(arg)
        elif opt in ("-h","--help"):
            usage()
            exit(0)    
        else :
            print ("\n Option not recognized %s\n\n" %(opt))
            usage()
            
           # assert False, "unhandled option"
            sys.exit(1)

    if argv_dict["bam_file"]== "" or argv_dict["output"] =="":
        print ("\n\tMissing required flags see usage\n")
        usage()
        exit(1)
    return argv_dict

if __name__ == "__main__":

    argv_dict = main_argv_parse(sys.argv)
    get_insert.parse_samfile(argv_dict["bam_file"], argv_dict["output"], argv_dict["depth"] )





