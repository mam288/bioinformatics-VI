"""
Solution to the Longest Shared Substring Problem.
Finding Mutations in DNA and Proteins (Bioinformatics VI) on Coursera.org
Week 1, code challenge #5
https://stepik.org/lesson/Suffix-Trees-296/step/6?unit=8994
"""

import networkx as nx
import wk1_04_longest_repeat as lr
import wk1_03_construct_suffix_trie as cst

def longest_shared_substring(input_data):
    '''
    Find the longest substring shared by two strings.
    
    Parameters
    --------
    input_data: input sequences (string)
    
    Return
    --------
    The longest substring that occurs in both input sequences (string)
    '''
    seqs = input_data.split()
    all_suffixes = [] 
    for seq in seqs:
        suffixes = cst.generate_suffixes(seq + '$')
        all_suffixes += suffixes
    trie = construct_suffix_trie_lss(all_suffixes)  # generate a trie for the combined suffixes of the input sequences
    common_substrings = repeats(trie)
    common_substrings = [x for x in common_substrings if all(x in seq for seq in seqs)]
    return max(common_substrings, key=len)  # return the longest string in common_substrings
    
def construct_suffix_trie_lss(patterns):
    '''
    Used with longest_shared_substring. Leaves are set to the next node number (not suffix starting index).
    ''' 
    trie = nx.DiGraph()
    trie.add_node(0)
    current_node = 0
    new_node_num = 1
    
    # create the suffix trie
    for pattern in patterns:
        sub_seq = ''
        for char in pattern:
            out_edges_list = trie.out_edges(current_node)
            out_edges_values = [trie.get_edge_data(incoming,outgoing)['val'] for (incoming,outgoing) in out_edges_list]
            if char == '$':
                trie.add_edge(current_node,new_node_num,{'val':'$'})  # set leaf value to the next node number
                new_node_num += 1
                current_node = 0
                break
            if char in out_edges_values:
                index = out_edges_values.index(char)
                current_node = out_edges_list[index][1]
            elif char not in out_edges_values:
                trie.add_edge(current_node,new_node_num,{'val':char})
                current_node = new_node_num 
                new_node_num += 1
    nodes = trie.nodes()
    
    # condense all non-branching portions of trie
    deg_1 =  [x for x in nodes if trie.in_degree(x) ==1 and trie.out_degree(x) == 1] 
    while deg_1 != []:
        for node in deg_1:
            pred = trie.predecessors(node)
            succ = trie.successors(node)
            sub_seq = trie[pred[0]][node]['val'] + trie[node][succ[0]]['val']
            trie.add_edge(pred[0],succ[0],{'val':sub_seq})
            trie.remove_edge(pred[0],node); trie.remove_edge(node,succ[0])
            trie.remove_node(node)
        nodes = trie.nodes()
        deg_1 =  [x for x in nodes if trie.in_degree(x) ==1 and trie.out_degree(x) == 1]
    return trie
    
def repeats(trie):
    '''Finds repeated substrings in a trie. Used with longest_shared_substring.'''
    edges = trie.edges()
    repeat_edge_list = []
    for (incoming,outgoing) in edges:
        if trie.out_degree(outgoing) > 1:
            path = lr.path_to_root(trie,outgoing)
            repeat_edge_list.append(path)
    return repeat_edge_list  

##############
if __name__ == "__main__":
    sample_input = '''TCGGTAGATTGCGCCCACTC
    AGGGGCTCGCAGTGTAAGAA'''
    print longest_shared_substring(sample_input)
