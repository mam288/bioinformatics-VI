"""
Solution to the Trie Construction Problem.
Finding Mutations in DNA and Proteins (Bioinformatics VI) on Coursera.org
Week 1, code challenge #1
https://stepik.org/lesson/Herding-Patterns-into-a-Trie-294/step/4?unit=8992
"""

import networkx as nx

def construct_trie(patterns):
    '''
    Solve the Trie Construction Problem.
    
    Parameters
    --------
    patterns: A collection of strings Patterns.
    
    Return
    ------
    Trie constructed from patterns
    
    Print output: The adjacency list corresponding to Trie(Patterns)
    '''
    patterns = patterns.split()
    patterns = [x + '$' for x in patterns]
    trie = nx.DiGraph()
    trie.add_node(0)
    current_node = 0
    new_node_num = 1
    
    # add each pattern to the trie
    for pattern in patterns:
        for char in pattern:
            out_edges_list = trie.out_edges(current_node)
            out_edges_values = [trie.get_edge_data(incoming,outgoing)['val'] for (incoming,outgoing) in out_edges_list]
            if char == '$':
                current_node = 0
                break
            if char in out_edges_values:
                index = out_edges_values.index(char)
                current_node = out_edges_list[index][1]
            elif char not in out_edges_values:
                trie.add_edge(current_node,new_node_num,{'val':char})
                current_node = new_node_num 
                new_node_num += 1
    
    # print the adjacency list for the trie
    for i in range(new_node_num):
        outgoing = trie.successors(i)
        for node in outgoing:
            char = trie.get_edge_data(i,node)['val']
            line1 = str(i) + '->' + str(node) + ':' + str(char)
            print line1
    return trie

###########################################
if __name__ == "__main__":
    patterns = '''ATAGA
    ATC
    GAT'''
    print construct_trie(patterns)
