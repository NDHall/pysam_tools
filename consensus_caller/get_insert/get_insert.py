



import pysam
from collections import Counter

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
            


                to_append="D%s" %( get_insert(pileupread,pileupcolumn.pos)[1:])
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
   

    for pileupcolumn in samfile.pileup():
        dict_str=''
        column_dict, contig = return_counter_for_col(pileupcolumn)
        ret_max, ret_support  = return_max_value(column_dict)
        print(contig, pileupcolumn.pos, ret_max, ret_support, column_dict)
        if ret_support >= min_req :
            if ret_max != "D":
                output_fasta_dict.setdefault(contig,["P"] * pileupcolumn.pos ).append(ret_max)
        elif ret_support < min_req :
            if ret_max != "D":
                output_fasta_dict.setdefault(contig,["P"] * pileupcolumn.pos ).append("N") 

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
 
            
    return output_fasta_dict


def parse_samfile(sam_file_name, out_file, min_req=20):

    samfile = pysam.AlignmentFile(sam_file_name, "rb" )

    output_fasta_dict = parse_pileup( samfile, min_req)
    f=open(out_file, 'w')
    
    for contig in output_fasta_dict:
        f.write(">%s\n" %(contig))
        f.write( "".join(output_fasta_dict[contig]) )
        f.write("\n")
    f.close()

            
    samfile.close()



if __name__  == "__main__"  :


    parse_samfile("/home/ndh0004/Documents/pysam_consensus/pysam_consensus/get_insert/equal_2.bam", "test", 20 )

