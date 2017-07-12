
def char_stats( char, charOccurrenceMap, totalCharCount ):
    """
    Explicity define, presence and absence of expected
    characters. In this case, I want to produce numbers
    for characters in set ["A","C","T","G","N","P"]. 
    I want to know when characters are not present and 
    to be able to print them to an expected. Column.
    """
    if char in charOccurrenceMap.keys():        
        charCount = charOccurrenceMap[char]
        relativeFrequency = charCount * 100.0 / totalCharCount
    else:
        charCount,relativeFrequency = [0,0]
    return [charCount, relativeFrequency ]



def printSequenceStatsCsv(fileName, sequenceName, charOccurrenceMap, totalCharCount):
    """
    Print details about a sequence to stdout.
    Called from inside parseFile().
    
    Keyword arguments:
        sequenceName: The name of the sequence to print
        charOccurrenceMap: A dict-like object that contains a char -->
                            number of occurrences mapping
        totalCharCount: The overall character count of the sequence
    """

    n_char, n_char_relative_freq = char_stats( "N", charOccurrenceMap, totalCharCount )
    p_char, p_char_relative_freq = char_stats( "P", charOccurrenceMap, totalCharCount )
    a_char, a_char_relative_freq = char_stats( "A", charOccurrenceMap, totalCharCount )
    t_char, t_char_relative_freq = char_stats( "T", charOccurrenceMap, totalCharCount )
    c_char, c_char_relative_freq = char_stats( "C", charOccurrenceMap, totalCharCount )
    g_char, g_char_relative_freq = char_stats( "G", charOccurrenceMap, totalCharCount ) 
    #print ("## file_name, seq_name, len, A, A_rel_freq, T, T_rel_freq, C, G_rel_freq, G, G_rel_freq, N, N_rel_freq, P, P_rel_freq" )
    print ( "%s,%s,%d,%d,%f%%,%d,%f%%,%d,%f%%,%d,%f%%,%d,%f%%,%d,%f%%" % ( fileName, sequenceName, totalCharCount, 
                                                                                                           a_char, a_char_relative_freq, 
                                                                                                           t_char, t_char_relative_freq, 
                                                                                                           c_char, c_char_relative_freq, 
                                                                                                           g_char, g_char_relative_freq,
                                                                                                           n_char, n_char_relative_freq,
                                                                                                           p_char, p_char_relative_freq))     
   
def columnNames():
    print ("file_name,seq_name,len,A,A_rel_freq,T,T_rel_freq,C,G_rel_freq,G,G_rel_freq,N,N_rel_freq,P,P_rel_freq" )
