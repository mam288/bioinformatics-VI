# -*- coding: utf-8 -*-
"""
Solution to the Inverse Burrows-Wheeler Transform Problem.
Finding Mutations in DNA and Proteins (Bioinformatics VI) on Coursera.org
Week 2, code challenge #3
https://stepik.org/lesson/The-First-Last-Property-and-Burrows-Wheeler-Inversion-299/step/10?unit=9001
"""

def bw_to_string(bw_seq):
    '''
    Reconstruct a string from its Burrows-Wheeler transform.
    
    Parameters
    --------
    bw_seq: burrows-wheeler sequence (string)
    
    Return
    --------
    original sequence used to create bw_seq (string)
    '''
    seq = ''
    first_column = ''.join(sorted(bw_seq))
    last_column = number_characters(bw_seq)
    first_column = number_characters(first_column)
    next_char = ''
    seq_index = 0
    while len(seq) < len(last_column):
        char,occurance = first_column[seq_index]
        bw_index = last_column.index((char,occurance))
        next_char = first_column[bw_index][0]
        seq += next_char
        seq_index = bw_index
    return seq
    
    
def number_characters(seq):
    '''
    Used with bw_to_string.
    Takes a given seq and returns a list of the same length with each element being a tuple (c,n). "c" is 
    the character at the matching index  in seq, and "n" is the occurance number of that character. 
    For example, ('A',3) would be the third occurance of 'A'.
    '''
    numbered_char_list = [('',-1) for i in range(len(seq))]
    completed_chars = []
    for char in seq:
        indices = [index for (index, character) in enumerate(seq) if character == char]
        if char not in completed_chars:
            count = 1
            for index in indices:
                numbered_char_list[index] = (char,count)
                count += 1
            completed_chars += [char]
    return numbered_char_list
      
###################################################################################
if __name__ == "__main__":
    sample_input = 'TTCCTAACG$A'
    challenge_data = 'TTGTCATCATGTGAACCTGCTTGGCAACTGACTATGGTGAGGCATAAGCCATCAACGAGCGCTAGCCCTATTCTGCATGTTGGACGTAAGTCCAGAAATCTCGTTAGGACGGAGCACCATCTCCAGACGGCTATTCGAAACCTTCGTTCCATGTCCGCAATAGCAAAAAAGTGTGCAAAATTATATTACGTCTAGGAATCTTCGACAGTGCCGCACCTTCTTTCGAAACCTGGCAGTAAGCCCTGCCCTACTAGGACGTTACTATT$GATCCCACATGCTACCTAACTCAGTTGCGGTACTGCCTCTGCTAGCTAAGGTGACTAAGAGGACTTATAAATTGAGAACAAGTCGATAACGTATGGGTAGGACCACCCCTCGAGCAGCATTTGCCGAGTCCTTTAAATGCACTGTCGACGGGTTCGAAATGCCACCTTGGTAGCGGCCAACCGCGACCTGGGGATGCTACATGAAAATTGGATAGGCGAGCATGCGAGAAATTGGAAACCCAGCTAAAATGAACCCAATACGGAACAACATGGACACTGCGGCACTGCGCTTACCTTTCCTTTCGTGGGCCAACACTTACCTAGAGACCACCCGCCTAAGCAGGTCTACATAAAGTCGAAGGTGGGCAGCCATATATGCCTACTCGACCTGGATATACTATGTCGAAACGCCGTTCCTTCGGTAAAACCTATAGCCCATGAGAAAGTATGGAATGGGCGTGTTGCAGATAGTCACGTATGATTTCAATTTGAAAGGCACCCAGGGCAGTTGCCTGTGCTTCGATGCCAATTTAGTGCGCAATCGCAAATGAGCAGCCAGACTTAGAATCATCGCCTAAGCGTCTGATACTCGCGTGGAAGGCCACAATGAGTGTAGCCACATGGAACC'
    print(bw_to_string(sample_input))
