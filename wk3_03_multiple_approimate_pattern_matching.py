# -*- coding: utf-8 -*-
"""
Solution to the Multiple Approximate Pattern Matching Problem.
Finding Mutations in DNA and Proteins (Bioinformatics VI) on Coursera.org
Week 3, code challenge #2
https://stepik.org/lesson/Epilogue-Mismatch-Tolerant-Read-Mapping-304/step/6?unit=9006
"""

import mutations as mu
import numpy as np
from numba import jit


def multiple_approximate_pattern_matching(seq,patterns,c,d):
    '''
    Multiple Approximate Pattern Matching Problem.
    
    Parameters
    --------
    seq: A DNA sequence (string)
    patterns: space delimited string of patterns (string)
    c: factor by which the checkpoint array is condensed (integer)
    d: maximum number of mismatches allowed (integer)
    
    Return
    --------
    All positions where one of the strings in Patterns appears as a substring of Text with at most d mismatches. 
    '''
    def create_count_dict(seq_array,c):
        '''
        Used with bw_better_matching. Creates count_dict which keeps track of the number of times a character appears in the array
        up to that point for each position in the array. Also creates first_occurances which keep strack of the index of the first
        occurance of each character in the array
        '''
        completed_chars = []
        first_occurances = {}
        array_len = (len(seq_array))/c + 1
        array_len += 1
        count_dict = {}
        for char in seq_array:
            if char not in count_dict:
                count_dict[char] = np.zeros((array_len,),dtype = int)
                indices = [index for (index, character) in enumerate(seq_array) if character == char]
                if indices == []:
                    continue
                else:
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

    lut = {' ':0,'A':1,  'C':2, 'G':3, 'T':4, '$':5}

    match_indices =  []
    bw_seq = mu.burrows_wheeler(seq + '$')
    bw_array = np.array([lut[x] for x in bw_seq]) # numpy array of bw_seq
    
    suffix_array = generate_suffix_array_list(seq + '$')
    
    
    first_column = ''.join(sorted(bw_seq))
    first_column_array = np.array([lut[x] for x in first_column]) # numpy array of firt_column
    
    count_dict_last,first_occurances_last_col = create_count_dict(bw_array,c)
    count_dict_first,first_occurances_first_col = create_count_dict(first_column_array,c)
    
    seq_array = np.array([lut[x] for x in seq]) # numpy array of seq
    
    patterns_array = np.array([lut[x] for x in patterns]) # numpy array of patterns
    patterns_zeros = np.where(patterns_array == 0)[0]
    patterns_zeros = np.concatenate([np.array([0]), patterns_zeros, np.array([len(patterns_zeros)])])
    def multiple_pattern_matching(seq,patterns,c,d,count_dict_last,first_occurances_first_col,suffix_array,bw_seq,first_column):
        '''
        Finds all starting positions in Text where a string from Patterns appears as a substring.
        '''
        
        def bw_matching_inner(first_occurances,count_dict,pattern,c,suffix_array):
            '''
            Returns the number of times pattern occurs in the sequence. 
            '''
            top = 0
            bottom = len(bw_seq) - 1
            while top <= bottom:
                if len(pattern) != 0:
                    current_char = int(pattern[-1])
                    pattern = pattern[:-1]
                    if current_char in bw_array[top:bottom+1]:
                        first_occurance_checkpoint = first_occurances[current_char]  # first occurance of char in the checkpoint array
                        first_occurance_seq = first_occurance_checkpoint*c + np.where(first_column[(first_occurance_checkpoint*c):(first_occurance_checkpoint*c+c)] == current_char)[0][0]        
                        top = first_occurance_seq + count_dict[current_char][top//c] + (bw_seq[(top//c)*c:top] == current_char).sum()
                        bottom = first_occurance_seq + count_dict[current_char][bottom//c] + (bw_seq[(bottom//c)*c:bottom+1] == current_char).sum() -1
                    else:
                        return -1
                else:                
                    pattern_match_locations = [suffix_array[i] for i in range(top,bottom+1)]
                    return pattern_match_locations
        all_match_locations = []

        for pattern in patterns:
            match_locations = bw_matching_inner(first_occurances_first_col,count_dict_last,str(pattern),c,suffix_array)
            if match_locations != -1:
                all_match_locations += match_locations
        return np.sort(all_match_locations)    
    
    
    # Use Just-in-time interpreter speedup from Numba with @jit decorator   
    @jit    
    def test_jit2(filtered_matches, seq_array, pattern, d):
        confirmed_matches = []
        for c_match_ind in xrange(len(filtered_matches)):
                match = filtered_matches[c_match_ind]
                match_seq = seq_array[match:match+len(pattern)]
                hamming_score = 0
                for i in xrange(len(match_seq)):
                    if match_seq[i] != pattern[i]:
                        hamming_score += 1
                if (hamming_score <= d):
                    confirmed_matches.append(match)
        return confirmed_matches
                    
    def divide_pattern(seq,kmer_len):
        substrings = []
        for i in range(len(seq)/kmer_len):
            i = i*kmer_len
            substring = seq[i:i+kmer_len]
            substrings.append(substring)
        if len(seq)%kmer_len != 0:
            np.concatenate([substrings[-1], seq[len(seq)/kmer_len*kmer_len:]])
        return substrings
    for i,e in enumerate(patterns_zeros):
        if i == len(patterns_zeros) - 1:
            break
        if i == len(patterns_zeros) - 2:
            pattern = patterns_array[patterns_zeros[i]+1:]
        elif i == 0:
            pattern = patterns_array[:patterns_zeros[i+1]]
        else:
            pattern = patterns_array[patterns_zeros[i]+1:patterns_zeros[i+1]]            
        kmer_len = len(pattern)/(d+1)
        seeds = divide_pattern(pattern,kmer_len)
        pattern_matches = []
        for i,seed in enumerate(seeds):
            potential_matches = multiple_pattern_matching(seq_array,seed,c,d,count_dict_last,first_occurances_first_col,suffix_array,bw_array,first_column_array)
            potential_matches = potential_matches - (kmer_len*i);
            potential_matches = potential_matches[potential_matches >= 0]

            
            
#            print potential_matches
            # get the seqeuence associated with each match
#            matches_seqs = [(seq_array[match:match+len(pattern)],match) for match in potential_matches if match+len(pattern) <= len(seq)]
            filtered_matches = potential_matches[potential_matches + len(pattern) <= len(seq)]
#            print pattern,matches_seqs[0][0]
            # remove any matches that have more than d mismatches
#            confirmed_matches = [n for (match_seq,n) in matches_seqs if [s[0]for s in difflib.ndiff(match_seq,seq)].count('+') <= d]
            
            confirmed_matches = test_jit2(filtered_matches, seq_array, pattern, d)
                                  
#            confirmed_matches = []
#            for c_match_ind in xrange(len(filtered_matches)):
#                match = filtered_matches[c_match_ind]
#                match_seq = seq_array[match:match+len(pattern)]
#                if (np.sum(match_seq != pattern) <= d):
#                    confirmed_matches.append(match)

#            confirmed_matches = [n for (match_seq,n) in matches_seqs if np.sum(match_seq != pattern) <= d] 
            
            pattern_matches += list(confirmed_matches) # confirmed matches for this seed
        match_indices += list(set(pattern_matches))
    match_indices.sort()
    match_list = ''
    for match in match_indices:
        match_list += str(match) + ' '
    write_file('results.txt',match_list)
    return match_list
def generate_suffixes_tuples(pattern):
    '''
    Generates all of the suffixes of pattern.
    INPUT: pattern (string)
    OUTPUT: list of all suffixes for pattern
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
    Generates a suffix array from a given sequence. The suffix array is a list of
    integers with each integer being the starting point of a suffix. The array is sorted 
    alphabetically by suffix.
    INPUT: seq(string)
    OUTPUT: index_list(list of integers)
    '''
    suffixes = sorted(generate_suffixes_tuples(seq),key=lambda x: x[0])    
    index_list = []
    for suffix in suffixes:
        index_list += [int(suffix[1])]
    return index_list        
        
    

def write_file(filename,content):
    with open(filename, 'a') as f:
        f.write(str(content))
#######################################################

#my_answer = (multiple_pattern_matching_mismatch(cd.challenge_data_seq,cd.challenge_data_patterns,100))
#my_answer = (multiple_pattern_matching(cd.extra_dataset_seq,cd.extra_dataset_patterns,5))
#print my_answer
#write_file('results_quiz.txt', my_answer)
#write_file('results_quiz.txt', '\n')
#write_file('results_quiz.txt', cd.a)
#correct_answer = cd.a.split()
#my_answer = my_answer.split()
sample_data_seq = 'ACATGCTACTTT'
sample_dataset_patterns = 'ATT GCC GCTA TATT'
#sample_dataset_patterns = 'ATT'
#start_time = time.time()
#print len(cd.extra_dataset_patterns.split())
#print multiple_approximate_pattern_matching(cd.challenge_data_seq,cd.challenge_data_patterns,1,2)
print (multiple_approximate_pattern_matching(sample_data_seq,sample_dataset_patterns,5,1))
#print("--- %s seconds ---" % (time.time() - start_time))
#print cd.extra_dataset_patterns[:57],len(cd.extra_dataset_patterns[:100])
#print cd.extra_dataset_seq[:1000],len(cd.extra_dataset_seq[:1000])
#print m.Approx_patterns(sample_data_seq,'ATT',1)
#for pattern in cd.extra_dataset_patterns.split():
#    
#write_file('results.txt', m.Approx_patterns(cd.extra_dataset_seq,pattern,cd.d))
#print len(my_answer),len(correct_answer)
#non_shared = set(my_answer)^set(correct_answer)

#print non_shared
#print len(non_shared)
#print t.divide_pattern_test('AAAAABBBBBCCCCCDDDDDEE',5)
#print generate_suffix_array_list('banana$')
#print mu.burrows_wheeler('CGTTTGCTAT$')
