"""
Implementation of the Viterbi algorithm solving the Decoding Problem.
Finding Mutations in DNA and Proteins (Bioinformatics VI) on Coursera.org
Week 4, code challenge #3
https://stepik.org/lesson/Hidden-Markov-Models-Code-Challenges-(Week-1)-11594/step/8?course=Stepic-Interactive-Text-for-Week-4&unit=9008
"""

import numpy as np
import networkx as nx
    
def viterbi(transition_matrix_input,emission_matrix_input,states,emission_letters,emission_string):
    '''
    Implement the Viterbi algorithm solving the Decoding Problem.
    
    Parameters
    ---------
    transition_matrix_input: initial Transition matrix (string)
    emission_matrix_input: initial Emission matrix (string)
    states:  Collection of states of the HMM (list of strings)      
    emission_letters: The HMM's alphabet Σ (list of strings)
    emission_string: A sequence of emitted symbols x = x1 . . . xn in an alphabet A generated by a k-state HMM with unknown transition and emission probabilities

    Return
    --------
     A path that maximizes the (unconditional) probability Pr(x, π) over all possible paths π. (string)
    '''
    def create_matrices(transition_matrix_input,emission_matrix_input,states,emission_letters):
        '''Convert the inputted transition and emission matrices to numpy arrays.'''
        num_emission_letters = len(emission_letters)
        num_states = len(states)
        transition_matrix = np.zeros((num_states,num_states))
        emission_matrix = np.zeros((num_states,num_emission_letters))
        
        # populate the transition matrix
        for i in range(1,num_states+1):
            row = transition_matrix_input.split('\n')[i]
            row = [float(x) for x in row.split()[1:]]
            transition_matrix[i-1] = row

        # populate the emission matrix
        for i in range(1,num_states+1):
            row = emission_matrix_input.split('\n')[i]
            row = [float(x) for x in row.split()[1:]]
            emission_matrix[i-1] = row
        return transition_matrix,emission_matrix
    
    def create_viterbi_graph(transition_matrix_input,emission_matrix_input,states,emission_letters,emission_string):
        '''Create a viterbi graph with weighted edges using the input parameters.'''
        viterbi_graph = nx.DiGraph()
        transition_matrix, emission_matrix = create_matrices(transition_matrix_input,emission_matrix_input,states,emission_letters)
        for i in range(len(emission_string)+1):
            if i == len(emission_string):
                index = 0
            else:
                index = emission_letters.index(emission_string[i])
            for outgoing_index in range(len(states)):
                outgoing_state =   states[outgoing_index]
                emission_weight = emission_matrix[outgoing_index][index]

                # calculate the weight for each edge (emission probability * transition probability) and populate the graph
                if i == 0:
                    incoming_state = ''
                    transition_weight = 0.5
                    viterbi_graph.add_edge(str(i),str(i+1) + outgoing_state,{'weight': emission_weight*transition_weight})
                else:
                    for incoming_index in range(len(states)):
                        incoming_state = states[incoming_index]
                        if i == len(emission_string):
                            transition_weight = emission_weight =  1
                            outgoing_state = ''
                            viterbi_graph.add_edge(str(i)  + incoming_state,str(i+1),{'weight': emission_weight*transition_weight})
                            continue 
                        transition_weight = transition_matrix[incoming_index][outgoing_index]                     
                        viterbi_graph.add_edge(str(i) + incoming_state,str(i+1) + outgoing_state,{'weight': emission_weight*transition_weight})
        return viterbi_graph
    
    def calculate_scores_viterbi(viterbi_graph):
        '''
        Calculate the score for each node in the viterbi graph
        '''
        viterbi_graph.node['0']['score'] = 1
        sorted_nodes = nx.topological_sort(viterbi_graph)
        backtrack = {}
        for node in sorted_nodes[1:]:
            predecessors = viterbi_graph.predecessors(node)
            high_score = 0
            backtrack_node = ''
            
            # calculate the path score coming from each predecessor and set the node score to the high score for all paths
            for predecessor in predecessors:
                score = viterbi_graph.node[predecessor]['score']*viterbi_graph[predecessor][node]['weight']
                if score > high_score:
                    high_score = score
                    backtrack_node = predecessor
            viterbi_graph.node[node]['score'] = high_score
            backtrack[node] = backtrack_node # document the prodecessor used to calculate the high score
        return backtrack
        
    def reconstruct_path(backtrack):
        '''
        Reconstruct the hidden path using the graph and backtrack matrix.
        '''
        last_node = str(len(emission_string) + 1)
        path = ''
        current_node = last_node
        predecessor = backtrack[current_node]
        while predecessor != '0':
            path = predecessor[-1] + path
            current_node = predecessor
            predecessor = backtrack[current_node]
        return path
    
    viterbi_graph = create_viterbi_graph(transition_matrix_input,emission_matrix_input,states,emission_letters,emission_string)
    backtrack = calculate_scores_viterbi(viterbi_graph)
    path = reconstruct_path(backtrack)
    return path
    
#######################################################################
if __name__ == "__main__":
    states = 'A B'.split()
    emission_string = 'xyxzzxyxyy'
    emission_letters = 'x y z'.split()
    transition_matrix = '''A	B
    A	0.641	0.359
    B	0.729	0.271'''
    emission_matrix = '''x	y	z
    A	0.117	0.691	0.192	
    B	0.097	0.42	0.483'''
    
    print viterbi(transition_matrix,emission_matrix,states,emission_letters,emission_string)
