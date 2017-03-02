"""
Solution to the Shortest Non-Shared Substring Problem.
Finding Mutations in DNA and Proteins (Bioinformatics VI) on Coursera.org
Week 1, code challenge #6
https://stepik.org/lesson/Suffix-Trees-296/step/6?unit=8994
"""

import wk1_05_longest_shared_substring as lss
import wk1_03_construct_suffix_trie as cst

def shortest_non_shared_substring(input_data):
    '''
    Shortest Non-Shared Substring Problem
    
    Parameters
    --------
    input_data: two amino acid sequences separated by a space (string)
        
    Return
    --------
    Shortest substring that is not shared between the input strings. (string)
    '''
    seqs = input_data.split()
    all_suffixes = []

    # generate a master list of the suffixes for all of the sequences
    for seq in seqs:  
        suffixes = cst.generate_suffixes(seq + '$')  # add a '$' to the end of the sequence and generate all suffixes
        all_suffixes += suffixes
        
    trie = lss.construct_suffix_trie_lss(all_suffixes)  # combine the suffixes from all input sequences and construct a trie
    possible_non_shared_substrings = possible_non_repeat_substrings(trie)
    non_shared_substrings = [x for x in possible_non_shared_substrings if not all(x in seq for seq in seqs)] # cull shared substrings
    min_len = len(min(non_shared_substrings, key = len))
    new_possible_non_shared_substrings = set([])
    
    # break each substring into substrings shorter than or equal to the length of the shortest substring in the original list
    for substring in non_shared_substrings:  
        new_substrings = generate_substrings(substring,min_len)  
        new_possible_non_shared_substrings = set(new_substrings)|set(new_possible_non_shared_substrings)
        
    new_possible_non_shared_substrings = list(new_possible_non_shared_substrings)
    new_non_shared_substrings = [x for x in new_possible_non_shared_substrings if not all(x in seq for seq in seqs)] # remove shared substrings
    return min(new_non_shared_substrings, key=len) 
    
    
def possible_non_repeat_substrings(trie):
    '''
    Used with shortest_non_shared_substring.
    Returns a set of substrings that might not be repeated in the input sequences. 
    '''
    nodes = trie.nodes()
    possible_non_repeat_substring_list = []

    # go through the nodes and find paths from each node to the prvious branching point
    for node in nodes:
        if trie.out_degree(node) == 0: # if node is a leaf
            path = path_to_prev_branch(trie,node) # append the path from the leaf to the previous branching point. 
            possible_non_repeat_substring_list.append(path)
            
    substring_list_no_dollarsign = [x for x in possible_non_repeat_substring_list if x[-1] != '$']  # list of substrings that do not end in $
    substring_list_remove_dollarsign = [x[:-1] for x in possible_non_repeat_substring_list if x[-1] == '$']  # list of substrings that previously ended in $, with the $ removed
    combined_substring_list = substring_list_no_dollarsign + substring_list_remove_dollarsign
    substring_list_no_blanks = [x for x in combined_substring_list if x != ''] # remove blank substrings
    substring_list_no_blanks = list(set(substring_list_no_blanks))  # remove repeats
    return substring_list_no_blanks
    
def path_to_prev_branch(trie,node):
    '''
    Used with shortest_non_shared_substring.
    Finds the path from the given node to the previous branching point in trie. 
    '''
    path = ''
    current_node = node
    while trie.out_degree(current_node) < 2 and current_node != 0:
        pred = trie.predecessors(current_node)
        edge_val = trie[pred[0]][current_node]['val']
        path = edge_val + path
        current_node = pred[0]
    return path

def generate_substrings(pattern,n):
    '''
    Generate all stubstrings up to length n for pattern.
    '''
    substrings = []
    for i in range(2,n+1):
        for j in range(len(pattern) - i +1):
            substring = pattern[j:j+i]
            substrings.append(substring)
    return substrings 

##############
if __name__ == "__main__":
    sample_input = '''CCAAGCTGCTAGAGG
    CATGCTGGGCTGGCT'''
    print shortest_non_shared_substring(sample_input)
