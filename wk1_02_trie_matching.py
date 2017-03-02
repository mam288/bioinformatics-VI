"""
Implementation of TrieMatching to solve the Multiple Pattern Matching Problem.
Finding Mutations in DNA and Proteins (Bioinformatics VI) on Coursera.org
Week 1, code challenge #2
https://stepik.org/lesson/Herding-Patterns-into-a-Trie-294/step/8?unit=8992
"""

import networkx as nx

def trie_matching(seq,patterns):
    '''
    Implement TrieMatching to solve the Multiple Pattern Matching Problem.
    Calculates starting positions in Text where a string from Patterns appears as a substring.
    
    Parameters
    --------
    seq: given DNA sequence (string)
    patterns: space delimited string of patterns (string)
    
    Return
    --------
    All starting positions in Text where a string from Patterns appears as a substring. (string with positions separated by spaces)
    '''
    trie = construct_trie(patterns)
    patterns = patterns.split()
    patterns = [x + '$' for x in patterns]
    seq += '$'
    current_index = 0
    start_index_list = ''
    while seq != '':
        pattern = prefix_trie_matching(seq,trie,patterns)
        if pattern != 'no matches found':
            start_index_list += str(current_index) + ' '
        seq = seq[1:]
        current_index += 1
    return start_index_list
    
def prefix_trie_matching(seq,trie,patterns):
    '''Traverses the trie to find locations where one of the given patterns is contained withing the given seq.'''
    seq_index = 0
    current_char = seq[seq_index]
    current_node = 0 # start at the root
    match = ''
    while True:    
        out_edges_list = trie.out_edges(current_node)
        out_edges_values = [trie.get_edge_data(incoming,outgoing)['val'] for (incoming,outgoing) in out_edges_list]                                
        if trie.in_degree(current_node) == 1 and trie.out_degree(current_node) == 0:  
            return match  # if you reach a leaf, return match
        elif current_char in out_edges_values: # if there is an outgoing edge equal to the next nucletide in the pattern
            match += current_char   # add the nucletide to 'match'
            
            #continue traversing the graph (update the index, char, and node)
            outgoing_edge_index = out_edges_values.index(current_char)
            current_node = out_edges_list[outgoing_edge_index][1] 
            seq_index += 1
            current_char = seq[seq_index] 
        else:
            return 'no matches found'
          
def construct_trie(patterns):
    '''Generates a trie from the given patterns.'''
    patterns = patterns.split()
    patterns = [x + '$' for x in patterns] # add a '$' to each pattern to indicate when we have reached the end of the pattern.
    trie = nx.DiGraph()
    trie.add_node(0)
    current_node = 0
    new_node_num = 1
    for pattern in patterns:
        for char in pattern:
            out_edges_list = trie.out_edges(current_node)
            out_edges_values = [trie.get_edge_data(incoming,outgoing)['val'] for (incoming,outgoing) in out_edges_list]
            if char == '$':  # if we reached the end of the pattern, set the current node back to 0 and move on to the next pattern
                current_node = 0
                break
            if char in out_edges_values:
                index = out_edges_values.index(char)
                current_node = out_edges_list[index][1]
            elif char not in out_edges_values:  # if the current node does not have an outgoing edge equal to the next nucletide in pattern
                trie.add_edge(current_node,new_node_num,{'val':char})  # add a new edge to the trie with the edge value equal to the next nucleotide in pattern
                current_node = new_node_num 
                new_node_num += 1
    for i in range(new_node_num):    # output the trie as a series of edges
        outgoing = trie.successors(i)
        for node in outgoing:
            char = trie.get_edge_data(i,node)['val']
            line1 = str(i) + '->' + str(node) + ':' + str(char)
            write_file('results_quiz.txt', line1)
    return trie


def write_file(filename,content):
    with open(filename, 'a') as f:
        f.write(content + '\n')    

###########################################
if __name__ == "__main__":
    seq_sample_input = 'AATCGGGTTCAATCGGGGT'
    patterns_sample_input = '''ATCG
    GGGT'''
    print trie_matching(seq_sample_input,patterns_sample_input)
