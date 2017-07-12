



import pysam
from collections import Counter
import logging

LEVELS = {'debug': logging.DEBUG,       
'info': logging.INFO,
'warning': logging.WARNING,
'error': logging.ERROR,
'critical': logging.CRITICAL}

logging.basicConfig(level=logging.WARNING)

def message(case, variables):
    if case=="cstring_error":
        message="""
        
        unrecognized cigar tuple value : %i
	accepted intergers: 4 = soft_clip
                            0 = match
			    1 = insertion
			    2 = deletion\n\n""" %(variables)
    elif case == "position_insert":
        read_pos, insert_ref_pos = variables
        message="""
        
        failure in cigar parsing addition!
        read_pos (%i) >insert_ref_pos (%i)
        an error has occured. Verify that
        all cigar tuples are correct, and 
        that no other weirdness has occured.

        \n""" %( read_pos, insert_ref_pos ) 


    elif case == "insert_seq_fail" :
        pos,ref,read = variables
        message="""
        No insertion sequence extracted for 
        at position %i in %s 
        for read : %s\n\n"""  %( pos, ref, read)



    elif case == "ref_skip":
        message="""
        A read with intron padding
        had been detected. This 
        function cannot parse reads
        with intron padding, Cigar
        string == "N": 
        read: %s \n\n""" %( variables)

    elif case == "complete" :
        bam, fasta, num_del, ins, N, P = variables 
        message="""
        Finished processing %s 
        Consensus sequence is written to %s:
        Total of %i deletions
        Total of %i insertions
        Total of %i bases below cutoff
        Total of %i bases with no coverage
        Breakdown by Contig.

""" %( bam, fasta, num_del, ins, N, P ) 

    elif case == "duplicate_contig" :
        message="""
        Duplicate contig %s discovered.\n\n""" %s( variables )



    return message


def get_insert( pileupread,col_number):
    insert_seq=None
    insert_ref_pos = col_number 
    insert_read_pos = pileupread.alignment.reference_start - 1
    read_coordinates=[0,0]
                
    for cstring,clen in  pileupread.alignment.cigartuples :
        #Accepted tuple values for cigar Identifiers
        soft_clip = 4
        match = 0
        insertion = 1
        deletion = 2
            

        # this conditional grabs the string when the postion in the 
        # cigar tuple list is equal to the position in the mpileup
        # and when the tuple is equal to an insertion.
        # 
        # the calculation of the start position is dropped by 1 so that it includes
        # the reference position. 

        if insert_read_pos == insert_ref_pos and cstring == insertion:
            read_coordinates[1]=read_coordinates[0]+clen
            start, stop = read_coordinates
            start=int(start -1) 
            stop=int(stop)
            insert_seq = pileupread.alignment.query_sequence[start : stop]
            break

        # this conditional calculates the correct starting coordinate on the 
        # read. It must include soft_clipping, matching and other insertions.
        # because it must allow other insertions preceding the position of the target
        # insertion. Thus conditional must be below, the conditional statment
        # that finally calculates stop coordinates. 


        if cstring in [soft_clip, match,  insertion]:
            read_coordinates[0] += clen

        # this conditional determines the relative position of the cigar tuple to the
        # reference sequence. When it is equal we are at the cigar tuple associated 
        # this mpileup column. Specifically, the cigar tuple that immediately proceeds
        # this mpileup column.
        if cstring in [0,2]:
            insert_read_pos += clen



        # tests, to ensure that only expected cigars are processed and that 
        # the cigar arthemtic is correct

        assert ( cstring in [soft_clip, match, deletion,insertion] ), message("cstring_error", cstring)
        assert(   insert_read_pos <= insert_ref_pos), message("position_insert",[insert_read_pos,insert_ref_pos] )  
    assert (insert_seq is not None), message("insert_seq_fail", [insert_ref_pos, pileupread.alignment.reference_name, pileupread.alignment.query_name ]) 
    return insert_seq               
    

def return_counter_for_col(pileupcolumn):

    # This function tallies, insertions, deletions, and mapped bases
    # deletions, are noted with D
    # insertions are returned as 5' match  +insertion_seq
    #  matches or snps returned as the base themselves. 
    contig=None
    summary_list=[]
    for pileupread in pileupcolumn.pileups:
        contig=pileupread.alignment.reference_name
        if not pileupread.is_del :
        # query position is None if is_del or is_refskip is set.
            if pileupread.indel >0 :

                # in this case indel seq must be called before single ref base
                # if not then this results in appending 2 possible seqs for same site. 

               summary_list.append(get_insert(pileupread,pileupcolumn.pos))
            else:
               summary_list.append( pileupread.alignment.query_sequence[pileupread.query_position])
            
        elif pileupread.is_del :
            if pileupread.indel >0 :
                
            # Currently handling deletions as they occur against the 
            # reference seqeunce. This is different from how insertions are 
            # handled. Insertions are handled as 1 block, meaning they are
            # are the result of 1 event. Deletion method assumes that deletions
            # will be consisent enough to provide the correct block like signal
            # even when each deletion against the reference is determined 
            # independently. If this creates severe systematic inconsistency.
            # it will have to be revisited. But if indel realignment is executed
            # such as in gatk, this should not be an issue.
            


                to_append="D%s" %( get_insert(pileupread,pileupcolumn.pos)[1:]) # this works because the raw coordinates return the base in the read before the deletion + the insertion. So that read has to be trimmed off.
            else :
                to_append="D"
            summary_list.append(to_append)
        assert( not  pileupread.is_refskip ), message("ref_skip", pileupread.query_name)
   
    return [( Counter(summary_list)), contig ]



def return_max_value(column_dict ):
    
    # python breaks ties by "randomness" of dict.
    # TODO: add a comparison to reference at same 
    # position and use the reference as a tie breaker.
    # additional work will have to be done if reference
    # is desired for tie breaker. 
    # this method is basically random given equal mapping
    # values. Looking at bowtie2 it does not always produce
    # identical mappings for identical, reads. While this is
    # is known, and often not problematic. It means that if 
    # additional steps like gatk realigner should be taken 
    # to maximize quality of alignment, since slight differences
    # can change the dominate read, using gatk will insure that
    # some degree of consistency is maintained across the whole
    # map.
    # example
    
    m=max(column_dict.values())

    ret_max =[k for k,v in column_dict.items() if v==m][0] 
    return [ ret_max, m ]



def parse_pileup( samfile, min_req=20):
    contig=None
    output_fasta_dict={}
    debug_pos=[]
    num_del = []
    num_ins = []
    prev_col = 0
    prev_contig = None
    for pileupcolumn in samfile.pileup():
        dict_str=''
        padding = ''
        column_dict, contig = return_counter_for_col(pileupcolumn)
        logging.debug("\n\nbegining of padding_loop",contig,  output_fasta_dict.keys())
        logging.debug("padding conditional", pileupcolumn.pos - prev_col >1 , prev_contig == contig , contig in output_fasta_dict.keys())
        if prev_contig != contig:
            logging.debug( "contig reset", prev_contig , contig)
            prev_contig = contig	
            padding = ''
            prev_col = pileupcolumn.pos
           

 


        elif  pileupcolumn.pos - prev_col >1 and prev_contig == contig and contig in output_fasta_dict.keys() :
            padding_int = pileupcolumn.pos - prev_col -1
            padding= "".join(["P"] * padding_int )
            logging.debug( "padding calculated", contig, padding, pileupcolumn.pos, prev_col )
 
        else:
            padding = ''


        # the above conditionals add internal padding.
        # briefly, 
	    # if checks that prev_col and mpileup.pos are
        #     on the same contig. If not there is no need for padding
        #     but a contig reset is required.
        # elif sets padding if step between positionsis greater than 1
        # else there is no need for padding.
        # padding is maintained as string and added to string character
        # returned by  function return_max_value().

           

        ret_max, ret_support  = return_max_value(column_dict)

        logging.debug("Begin append conditionals",contig, pileupcolumn.pos, ret_max, ret_support, column_dict)
        logging.debug( "Contig",contig, "padding", padding)
        if ret_support >= min_req :
            if "D" not in ret_max:
                output_fasta_dict.setdefault(contig,["P"] * pileupcolumn.pos ).append( padding + ret_max)
                #print ( "if", contig, ["P"] * pileupcolumn.pos, pileupcolumn.pos )
                logging.debug("No-Delete", contig, pileupcolumn.pos, min_req,"min_req")
                if len(ret_max)> 1:
                    num_ins.append(contig)
            elif "D" in ret_max and len(ret_max) >1 :
                num_del.append(contig) 
                logging.debug( "deletion", contig, ["P"] * pileupcolumn.pos, pileupcolumn.pos )
                output_fasta_dict.setdefault(contig,["P"] * pileupcolumn.pos ).append( padding + ret_max[1:])
                logging.debug("Delete", contig, pileupcolumn.pos,  min_req,"min_req")
            elif "D" in ret_max and len(ret_max)==1:
                num_del.append(contig)
        elif ret_support < min_req :
            output_fasta_dict.setdefault(contig,["P"] * pileupcolumn.pos ).append( padding + "N") 
            #print("Low-support", contig, pileupcolumn.pos, output_fasta_dict[contig])
            logging.debug("Low-support", contig, prev_contig, pileupcolumn.pos, prev_col,  min_req,"min_req",padding,
                  ( pileupcolumn.pos - prev_col >1) ,(  prev_contig == contig )
                          )

                    # even if it was an insertion, this 
                    # okay, because the insertion is not 
                    # supported and reference pos is used 
                    # instead.
                    #
                    # Also, note that the first instance of the 
                    # contig region is set to pad the alignment
                    # fully mapped alignment will start at pos =0
                    # if there is an unmapped region that preceeds
                    # the first position it is padded with character 
                    # P. This is to distguish it from the N used for
                    # bases that do not meet the minimum requirement.
        prev_col = pileupcolumn.pos
        logging.debug("Keys at for output_fasta_dict",output_fasta_dict.keys(),"\nreturned for contig",contig)
 
            
    return [output_fasta_dict , Counter(num_del), Counter(num_ins)]


def parse_samfile(sam_file_name, out_file, min_req=20):

    samfile = pysam.AlignmentFile(sam_file_name, "rb" )

    output_fasta_dict ,num_del, num_ins = parse_pileup( samfile, min_req)
    f=open(out_file, 'w')
    P = 0
    N = 0
    padding_ins_dict={}
    for contig in output_fasta_dict:
        local_p = 0 
        local_n = 0
        contig_len = 0 
        f.write(">%s\n" %(contig))
        f.write( "".join(output_fasta_dict[contig]) )
        f.write("\n")
        local_p = "".join(output_fasta_dict[contig]).count("P")
        local_n = "".join(output_fasta_dict[contig]).count("N")
        contig_len = len( output_fasta_dict[contig] )
        assert (contig not in padding_ins_dict), message("duplicate_contig", contig)
        padding_ins_dict[contig] = [local_p, local_n, contig_len ]
        P += local_p
        N += local_n 
    f.close()

            
    samfile.close()

    print( message("complete", [ sam_file_name, out_file, len(num_del.values()), len(num_ins.values()), N, P]))
    print ("\n\n\t\tcontig\tlen\tins\tdel\tlow-cov\tno-cov" )
    for contig in padding_ins_dict :
       P, N, contig_len = padding_ins_dict[contig]

       if contig in num_ins.keys():
           num_ins_ret = num_ins[contig]
       else:
           num_ins_ret = 0 
       if contig in num_del.keys():
           num_del_ret = num_del[contig]
       else :
           num_del_ret = 0    

       print ( "\t\t%s\t%i\t%i\t%i\t%i\t%i" %( contig, contig_len, num_ins_ret, num_del_ret, N, P ) )

 
    


if __name__  == "__main__"  :


    parse_samfile("/home/ndh0004/Documents/pysam_consensus/pysam_consensus/get_insert/equal_2.bam", "test", 20 )

