"""
Solution to the Suffix Tree Construction Problem.
Finding Mutations in DNA and Proteins (Bioinformatics VI) on Coursera.org
Week 1, code challenge #3
https://stepik.org/lesson/Suffix-Trees-296/step/4?unit=8994
"""

import networkx as nx

def construct_suffix_trie(seq):
    '''
    Suffix Tree Construction Problem.
    Constructs a trie from all of the suffixes of a given sequence. Condense all non-branching stretches of the trie.
    
    Parameters
    --------
    seq: DNA sequence (string) 
    
    Return
    --------
    Trie (NetworkX DiGraph)
    
    Print Output: The edge labels of SuffixTree(Text). 
    '''
    patterns = generate_suffixes(seq)
    trie = nx.DiGraph()  
    trie.add_node(0)
    current_node = 0
    new_node_num = 1
    
    # create a suffix trie
    for pattern in patterns:
        sub_seq = ''
        for char in pattern:
            out_edges_list = trie.out_edges(current_node)
            out_edges_values = [trie.get_edge_data(incoming,outgoing)['val'] for (incoming,outgoing) in out_edges_list]
            if char == '$':
                trie.add_edge(current_node,str(len(seq) - len(pattern)),{'val':'$'})
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
 
    # generate a list of all nodes along nodes along non-branching segments of trie (in-degree and out-degree of 1)
    deg_1 =  [x for x in nodes if trie.in_degree(x) ==1 and trie.out_degree(x) == 1]  
              
    # condense all non-branching portions of trie
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
    edges = trie.edges()
    for (incoming,outgoing) in edges:  # print all of the edge labels for trie
#        write_file('results_quiz.txt', trie[incoming][outgoing]['val'])
         print (trie[incoming][outgoing]['val'])
    return trie

def generate_suffixes(pattern):
    '''
    Generates all of the suffixes of pattern.
    '''
    suffixes = []
    while pattern != '':
        suffixes.append(pattern)
        pattern = pattern[1:]
    return suffixes  
    
##############
if __name__ == "__main__":
    sample_input = 'ATAAATG$'
    construct_suffix_trie(sample_input)
