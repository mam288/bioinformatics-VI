# -*- coding: utf-8 -*-
"""
Solution to the Multiple Pattern Matching Problem.
Finding Mutations in DNA and Proteins (Bioinformatics VI) on Coursera.org
Week 3, code challenge #2
https://stepik.org/lesson/Burrows-and-Wheeler-Set-Up-Checkpoints-303/step/4?unit=9005
"""
import mutations as mu
import numpy as np

def multiple_pattern_matching(seq,patterns,c):
    '''
    Multiple Pattern Matching Problem. Calculates the number of times each pattern occurs 
    in the original sequence for bw_seq. Faster and less memory than the original bw_matching.
    
    Parameters
    --------
    seq: A DNA sequence (string)
    patterns: space delimited string of patterns (string)
    c: factor by which the checkpoint array is condensed (integer)
        
    Return
    --------
    All starting positions in Text where a string from Patterns appears as a substring. (space separated string)
    '''
        
    def bw_matching_inner(first_occurances,count_dict,pattern,c,suffix_array):
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
                    first_occurance_checkpoint = first_occurances[current_char]  # first occurance of char in the checkpoint array
                    first_occurance_seq = first_occurance_checkpoint*c + first_column[(first_occurance_checkpoint*c):(first_occurance_checkpoint*c+c)].index(current_char)           
                    top =  first_occurance_seq + count_dict[current_char][top//c] + bw_seq[(top//c)*c:top].count(current_char)
                    bottom = first_occurance_seq + count_dict[current_char][bottom//c] + bw_seq[(bottom//c)*c:bottom+1].count(current_char) -1
                else:
                    return -1
            else:                
                pattern_match_locations = [suffix_array[i] for i in range(top,bottom+1)]
                return pattern_match_locations
                
    bw_seq = mu.burrows_wheeler(seq + '$')
    count_dict_last,first_occurances_last_col = create_count_dict(bw_seq,c)
    first_column = ''.join(sorted(bw_seq))
    count_dict_first,first_occurances_first_col = create_count_dict(first_column,c)
    suffix_array = generate_suffix_array_list(seq + '$')
    patterns = patterns.split()
    all_match_locations = []

    # find the match locations for each pattern and add them to the master list of locations
    for pattern in patterns:
        match_locations = bw_matching_inner(first_occurances_first_col,count_dict_last,pattern,c,suffix_array)
        if match_locations != -1:
            all_match_locations += match_locations
            
    sorted_all = sorted(all_match_locations)
    match_locations_string = ''
    for location in sorted_all:
        match_locations_string += str(location) + ' '
    return match_locations_string

    
def generate_suffixes_tuples(pattern):
    '''
    Generates all of the suffixes of pattern.
    '''
    suffixes = []
    index = 0
    while pattern != '':
        suffixes.append((pattern,index))
        index += 1
        pattern = pattern[1:]
    return suffixes
    
def generate_suffix_array_list(seq):
    '''
    Generates a suffix array from a given sequence. The suffix array is a string of
    integers with each integer being the starting point of a suffix. The array is sorted 
    alphabetically by suffix.
    '''
    suffixes = sorted(generate_suffixes_tuples(seq),key=lambda x: x[0])
    index_list = []
    for suffix in suffixes:
        index_list += [int(suffix[1])]
    return index_list

def construct_suffix_array(word):
    '''Constructs a suffix array from the given word.'''
    word += ['', '$'][word[-1] != '$'] # Check that the word ends in '$'.
    suffix_comp = lambda i,j: [1, -1][word[i] < word[j]] if word[i] != word[j] else suffix_comp(i+1,j+1)
    suffix_array = sorted(xrange(len(word)), cmp=suffix_comp)
    return suffix_array

def create_count_dict(seq,c):
    '''
    Used with bw_better_matching. Creates count_dict which keeps track of the number of times a character appears in the array
    up to that point for each position in the array. Also creates first_occurances which keep strack of the index of the first
    occurance of each character in the array
    '''
    completed_chars = []
    first_occurances = {}
    array_len = (len(seq))/c + 1
    array_len += 1
    count_dict = {'A':np.zeros((array_len,),dtype = int),'C':np.zeros((array_len,),dtype = int),'G':np.zeros((array_len,),dtype = int),'T':np.zeros((array_len,),dtype = int),'$':np.zeros((array_len,),dtype = int)}
    for char in seq:
        indices = [index for (index, character) in enumerate(seq) if character == char]
        first_occurances[char] = indices[0]/c
        if char not in completed_chars:
            count = 0
            for i in range(len(indices)):
                index = indices[i]
                if i >0:
                    count_dict[char][(indices[i-1])/c +1:(indices[i])/c+1] = count
                count += 1
            count_dict[char][(indices[i])/c+1:] = count
            completed_chars += [char]

    return count_dict, first_occurances
    

def write_file(filename,content):
    with open(filename, 'a') as f:
        f.write(str(content))
#######################################################


sample_data_seq = 'AATCGGGTTCAATCGGGGT'
sample_data_patterns= '''ATCG
GGGT'''
print (multiple_pattern_matching(sample_data_seq,sample_data_patterns,5))
