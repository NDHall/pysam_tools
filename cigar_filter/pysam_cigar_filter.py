import pysam


"""
    This script is for filtering out soft-clipped reads that are speciously mapped.
    It is of particular concern especially if no attempt has been made to filter 
    reads during mapping for duplicate regions such as are found between chloroplast
    and mitochondrial regions. Even if filtering is applied, this filter will remove
    soft clipped reads that are suspcious.

"""
import logging
import getopt
import sys



LEVELS = {'debug': logging.DEBUG,
          'info': logging.INFO,
          'warning': logging.WARNING,
          'error': logging.ERROR,
          'critical': logging.CRITICAL}

logging.basicConfig(level=logging.WARNING)

logging.debug("\tLogging level is set to debug\n")

def usage():
    message="""
    pysam_cigar_filter.py  -s or --sam_file <sam or bam file> 
                           -l --length <Minimum acceptable match length int default=20> 
                           -o --out_prefix <out_prefix> 
                           -h --help [returns help option]
"""
    print (message)
    return message

def messages(case,variable):
    if case == "Reference_Trouble":
        message="""

    There are repeated reference names in sam header.
    This must be fixed. The repeated name is : %s

""" %( variable )

    elif case == "mapped_reads":
        total_read_counter , mapped_read_counter = variable
        message="""
        
        
        Total of %i reads processed
        Total of %i reads excluded because of soft clipping
        Total of %i reads included.
        
        """ %( total_read_counter, (total_read_counter - mapped_read_counter), mapped_read_counter)

    print (message)

def main_argv_parse(argv):
    logging.debug("inside main_argv_parse sys.argv: %s" %(argv[1:]))
    argv_dict={"output_prefix":'' ,"sam_file":'',"length":20};
    if len(argv[1:]) == 1 and ( argv[1] == '-h' or argv[1]== '--help'):
        usage()
        exit(0)
        
    try:
        opts, args = getopt.gnu_getopt(argv[1:],"h:s:l:o:",["help=","output_prefix=","sam_file=","length="])
        logging.debug("opts: %s\nargs: %s "%(opts, args))
    except getopt.error:
        logging.critical("getopt.error argv: %s" %(" ".join(argv)))
        usage()
        sys.exit(2)
    for opt , arg in opts:
        if opt in ("-o","--output_prefix"):
            argv_dict["output_prefix"]=arg
        elif opt in ("-s","--sam_file"):
            argv_dict["sam_file"]=arg
        elif opt in ("-l","--length"):
            argv_dict["length"]=int(arg)
        elif opt in ("-h","--help"):
            usage()
            exit(0)    
        else :
            print ("\n Option not recognized %s\n\n" %(opt))
            usage()
            
           # assert False, "unhandled option"
            sys.exit(1)

    return argv_dict



def alignment_file(sam):
    ref_len_dict={}
    # max_ref_pos={"dummy_name":["current_max_pos",["list_of_reads..."]]}
    # min_ref_pos={"dummy_name":["current_max_pos",["list_of_reads..."]]}
    max_ref_pos={}
    min_ref_pos={}   
    if sam.split(".")[-1] =="sam" :
        sam_file=pysam.AlignmentFile(sam, 'r')
    elif sam.split(".")[-1] =="bam" :
        sam_file=pysam.AlignmentFile(sam, 'rb' )
    else:
        print("\n\n\tNo sam or bam extension detected for %s\n\n" %(sam))
        usage()
        exit(1)
    for SQ in sam_file.header["SQ"] :
        #SQ["SN"] is the reference name in the header dict.
        #SQ["LN"] is the length.
        #Make a dictionary of the expected last position for reads to map to.
        assert (SQ["SN"] not in ref_len_dict) , messages("Reference_Trouble", SQ["SN"]) 
        ref_len_dict[SQ["SN"]] = SQ["LN"]
        max_ref_pos[SQ["SN"]] = [20,[]]
        min_ref_pos[SQ["SN"]]=[1000,[]]
        #since max_ref_pos is just a set of names we can use it for
        # for both max_ref_pos and min_ref_pos dictionaries.

    return [sam_file, ref_len_dict, max_ref_pos, min_ref_pos ]

def cigar_read_filter(read, length,ref_len_dict , max_ref_pos, min_ref_pos, cutoff=5):
    """
    bam tuple ids used by pysam.
    right now we are using tupples 
     with ID 4 or 0
    cigartuples returns list of tuples with 
     Formatted (ID_Field, Len)
        M    BAM_CMATCH    0
        I    BAM_CINS    1
        D    BAM_CDEL    2
        N    BAM_CREF_SKIP    3
        S    BAM_CSOFT_CLIP    4
        H    BAM_CHARD_CLIP    5
        P    BAM_CPAD    6
        =    BAM_CEQUAL    7
        X    BAM_CDIFF    8
        B    BAM_CBACK    9

    check that read is greater than or equal to
     minimum allowable length.\
    """
    ret_read=False
    if length <= read.reference_length :
        soft_clip=4
        match=0
        cigar=read.cigartuples
        start=cigar[0]
        end=cigar[-1]
        # soft clippling by definition occurs at either end.
        # if it occurs at both ends it will be excluded as poor quality mapping.
        # soft clipping will only be allowed on the first or last mapped base pair. 
        if (start[0] == soft_clip and end[0]== soft_clip ):
            pass
                # if soft clipping is at the end for the last base, we want to allow allow that for last base. 
        elif end[0] == soft_clip :
            if end[1] < cutoff:
                ret_read=True
            elif max_ref_pos[read.reference_name][0]< read.reference_end:
                 max_ref_pos[read.reference_name]=[read.reference_end,[read]]
            elif max_ref_pos[read.reference_name][0]== read.reference_end:
                max_ref_pos[read.reference_name][1].append(read)
            else:
                pass
            # if the read is soft clipped at the 3prime end less than the
            # the maximum there is nothing to do with.
        
        # if reads do not start mapping at the first base of the contig we want to be able to catch
        # the reads that are the very first ones mapped and allow soft clipping.
       
        elif start[0]== soft_clip :
            if start[1] < cutoff:
                ret_read=True
            elif  read.reference_name in min_ref_pos:
                if read.reference_start == 0:
                    ret_read=True
 
                    del  min_ref_pos[read.reference_name]
                    logging.debug("Absolute Minimum Found for %s == 0" %(read.reference_name))

                elif read.reference_start < min_ref_pos[read.reference_name][0]:

                    logging.debug("New Minimum found for %s == new %i old %i " %(read.reference_name, read.reference_start, min_ref_pos[read.reference_name][0]))
                    min_ref_pos[read.reference_name]=[read.reference_start,[read]]
                elif  read.reference_start == min_ref_pos[read.reference_name][0] :
                   # print(read)
                    min_ref_pos[read.reference_name][1].append(read)
            elif read.reference_name not in min_ref_pos and read.reference_start == 0 : 
                ret_read=True
                
        else:
            ret_read=True

        if ret_read is True:
       
            return read 
        else :

            return None


def out_from_read_dict(out, read_dict, mapped_read_counter):
    for contigs in read_dict:
        for reads in read_dict[contigs][1]:
            mapped_read_counter+=1
            out.write(reads)
    return mapped_read_counter

    
            



def soft_clip_filter(sam_file, out, length, ref_len_dict, max_ref_pos, min_ref_pos ): 
    """ 
        This function iterates through the sam file and returns reads that are soft clipped on either end of
        the alignment. It also allows for one end to be soft-clipped up to 5 bp, but not both ends. This is 
        to allow for the case where not all of an adapter was removed. If a global alignment is forced in 
        this situation, then it introduces artificial noise into the alignment resulting in poorly called 
        bases.
    """ 
    mapped_read_counter=0
    total_read_counter=0

    for read in sam_file:
        total_read_counter+=1
        out_read =cigar_read_filter(read, length, ref_len_dict, max_ref_pos, min_ref_pos )
        if out_read is not None != 0 :
            logging.debug( "read passed")
            # logging.debug(read.cigartuples, read.reference_start, read.reference_end)
            out.write(out_read)
            mapped_read_counter+=1


    mapped_read_counter = out_from_read_dict(out, min_ref_pos, mapped_read_counter)
    mapped_read_counter = out_from_read_dict(out, max_ref_pos,mapped_read_counter)

    messages("mapped_reads",[total_read_counter,mapped_read_counter])
   

        
    
            
         



if __name__ == "__main__" :
    argv=main_argv_parse(sys.argv)

    logging.debug("argv :")
    logging.debug(argv)


    sam_file , ref_len_dict, max_ref_pos, min_ref_pos = alignment_file(argv["sam_file"])
    logging.debug("argv[\"length\"]")
    logging.debug(argv["length"])

    out_bam=pysam.AlignmentFile(str(argv["output_prefix"])+".bam", "wb", header=sam_file.header)
    
    soft_clip_filter( sam_file, out_bam ,argv["length"], ref_len_dict, max_ref_pos, min_ref_pos)
      
 
    out_bam.close()
    sam_file.close()
    

    pysam.sort("-o", "sorted_"+str(argv["output_prefix"])+".bam" , str(argv["output_prefix"])+".bam" )


    
