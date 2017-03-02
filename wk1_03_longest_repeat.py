# -*- coding: utf-8 -*-
"""
Solution to the Longest Repeat Problem.
Finding Mutations in DNA and Proteins (Bioinformatics VI) on Coursera.org
Week 1, code challenge #4
https://stepik.org/lesson/Suffix-Trees-296/step/5?unit=8994
"""

import wk1_03_construct_suffix_trie as cst

def longest_repeat(input_string):
    '''
    Find the longest repeat in a string.
    
    Parameters
    --------
    input_string: A DNA sequence (string)
    
    Return
    --------
    A longest substring of Text that appears in Text more than once. (string)
    '''
    trie = cst.construct_suffix_trie(input_string)
    edges = trie.edges()
    repeat_substrings = []   # list of substrings that are repeated in the input_string
    for (incoming,outgoing) in edges:  
        if trie.out_degree(outgoing) > 1: # if outgoing node branches (indicating the substring before the branch is repeated)
            path = path_to_root(trie,outgoing)  # find the path from the outgoing node to the root
            repeat_substrings.append(path)    # add the path to the list of repeat edges
    longest_repeat = max(repeat_substrings, key=len)  
    return longest_repeat

    
def path_to_root(trie,node):
    '''
    Finds the path from the given node to the root of trie. Path is outputted as string consisting of all
    the edge values along the path added together. 
    '''
    path = ''
    current_node = node
    while current_node != 0:
        pred = trie.predecessors(current_node)
        edge_val = trie[pred[0]][current_node]['val']
        path = edge_val + path
        current_node = pred[0]
    return path  

##############

sample_input = 'ATATCGTTTTATCGTT$'
print(longest_repeat(sample_input))
