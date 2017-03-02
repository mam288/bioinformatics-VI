"""
Implementation of BetterBWMatching.
Finding Mutations in DNA and Proteins (Bioinformatics VI) on Coursera.org
Week 3, code challenge #1
https://stepik.org/lesson/Speeding-Up-Burrows-Wheeler-Pattern-Matching-301/step/1?unit=9003
"""

import numpy as np

def better_bw_matching(bw_seq,patterns):
    '''
    Implement BetterBWMatching.
    Faster and less memory than the original bw_matching.
     
    Parameters
    --------
    bw_seq: Burrows-Wheeler transformation of A DNA sequence (string)
    patterns: space delimited string of patterns (string)
    
    Return
    A list of integers, where the i-th integer corresponds to the number of substring matches 
    of the i-th member of Patterns in Text. (string with integers separated by spaces)
    '''
        
    def bw_matching_inner(first_occurances,count_dict,pattern):
        '''
        Returns the number of times pattern occurs in the sequence. 
        '''
        top = 0
        bottom = len(bw_seq) - 1
        while top <= bottom:
            if pattern != '':
                current_char = pattern[-1]
                pattern = pattern[:-1]

                # calculate the new top and bottom rows for the range
                if current_char in bw_seq[top:bottom+1]:
                    top = first_occurances[current_char] + count_dict[current_char][top]
                    bottom = first_occurances[current_char] + count_dict[current_char][bottom + 1] -1
                else:
                    return 0
            else:
                return bottom - top + 1
                
    count_dict_last,first_occurances_last_col = create_count_dict(bw_seq)
    first_column = ''.join(sorted(bw_seq))
    count_dict_first,first_occurances_first_col = create_count_dict(first_column)
    patterns = patterns.split()
    occurances_string = ''
    for pattern in patterns:
        num_occurances = bw_matching_inner(first_occurances_first_col,count_dict_last,pattern)
        occurances_string += str(num_occurances) + ' '
    return occurances_string
    
    
def create_count_dict(seq):
    '''
    Used with bw_better_matching. Creates count_dict which keeps track of the number of times a character appears in the array
    up to that point for each position in the array. Also creates first_occurances which keep strack of the index of the first
    occurance of each character in the array.
    '''
    completed_chars = []
    first_occurances = {}
    count_dict = {'A':np.zeros((len(seq)+1,),dtype = int),'C':np.zeros((len(seq)+1,),dtype = int),'G':np.zeros((len(seq)+1,),dtype = int),'T':np.zeros((len(seq)+1,),dtype = int),'$':np.zeros((len(seq)+1,),dtype = int)}
    for char in seq:
        indices = [index for (index, character) in enumerate(seq) if character == char]
        first_occurances[char] = indices[0]
        if char not in completed_chars:
            count = 0
            for i in range(len(indices)):
                index = indices[i]
                if i >0:
                    count_dict[char][indices[i-1]+1:indices[i]+1] = count
                count += 1
            count_dict[char][indices[i]+1:] = count
            completed_chars += [char]
    return count_dict, first_occurances
    


#####################################################################################################

if __name__ == "__main__":
    sample_input_bw = 'GGCGCCGC$TAGTCACACACGCCGTA'
    sample_input_patterns = 'ACC CCG CAG'
    print (better_bw_matching(sample_input_bw,sample_input_patterns))
