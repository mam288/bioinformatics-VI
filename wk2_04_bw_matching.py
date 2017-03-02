"""
Implementation of BWMatching.
Finding Mutations in DNA and Proteins (Bioinformatics VI) on Coursera.org
Week 2, code challenge #4
https://stepik.org/lesson/Pattern-Matching-with-the-Burrows-Wheeler-Transform-300/step/8/toc?unit=9002
"""

def bw_matching(bw_seq,patterns):
    '''
    Implement BWMatching.
    Calculates the number of times each pattern occurs in the original sequence for bw_seq. 
    
    Parameters
    --------
    bw_seq: given bw_seq (string)
    patterns: space delimited string of patterns (string)
    
    Return
    --------
    A list of integers, where the i-th integer corresponds to the number of substring matches of the i-th member of Patterns
    in Text. (string with integers separated by spaces)
    '''
    def bw_matching_inner(first_column,last_column,pattern,last_to_first):
        '''Returns the number of times pattern occurs in the sequence.'''
        top = 0
        bottom = len(bw_seq) - 1
        while top <= bottom:
            if pattern != '':
                current_char = pattern[-1]
                pattern = pattern[:-1]

                # find the indices of the entries in last_column that match the current_char and are within the set top and bottom limits
                matching_indices = [index for (index,entry) in enumerate(last_column) if last_column[index][0] == current_char and index >= top and index <= bottom]
                if matching_indices == []:
                    return 0
                else:
                    top_index = matching_indices[0]
                    bottom_index = matching_indices[-1]
                    top = last_to_first[top_index]
                    bottom = last_to_first[bottom_index]
            else:
                return bottom - top + 1
    last_column = number_characters(bw_seq)
    first_column = ''.join(sorted(bw_seq))
    first_column = number_characters(first_column)
    last_to_first = [first_column.index(e) for e in last_column]
    patterns = patterns.split()
    occurances_string = ''
    
    # for each pattern - find the occurances in the given string and add to the master list
    for pattern in patterns:
        num_occurances = bw_matching_inner(first_column,last_column,pattern,last_to_first)
        occurances_string += str(num_occurances) + ' '
        
    return occurances_string
    
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
    
#####################################################################################################
if __name__ == "__main__":
    sample_input_bw = 'TCCTCTATGAGATCCTATTCTATGAAACCTTCA$GACCAAAATTCTCCGGC'
    sample_input_patterns = 'CCT CAC GAG CAG ATC'
    print (bw_matching(sample_input_bw,sample_input_patterns))