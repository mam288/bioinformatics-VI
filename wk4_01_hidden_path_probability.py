# -*- coding: utf-8 -*-
"""
Solution to the Probability of a Hidden Path Problem.
Finding Mutations in DNA and Proteins (Bioinformatics VI) on Coursera.org
Week 4, code challenge #1
https://stepik.org/lesson/Hidden-Markov-Models-Code-Challenges-(Week-1)-11594/step/2?course=Stepic-Interactive-Text-for-Week-4&unit=9008
"""

import numpy as np
import networkx as nx


    
def hidden_path_probability(input_matrix,states,path):
    '''
    Probability of a Hidden Path Problem.
    
    Parameters
    --------
    input_matrix: Emission matrix (string)
    states:  Collection of states of the HMM (list of strings)     
    path: hidden path π
   
    Return
    --------
    The probability of this path, Pr(π). (floating point number)
    '''
    transition_dict = create_dict_hmm(input_matrix,states)
    product = 0.5
    for i in range(1,len(path)):
        outgoing_index = states.index(path[i])
        weight = transition_dict[path[i-1]][outgoing_index]
        product = product*weight
    return product

def create_dict_hmm(input_matrix,states):
    '''
    Convert the input matrix into a dicionary.
    '''
    num_states = len(states)
    input_matrix_array = np.array(input_matrix.split())
    transition_matrix = np.array([float(x) for x in input_matrix_array if x not in states])
    transition_dict = {}
    for i,state in enumerate(states):
        values = transition_matrix[i*num_states:i*num_states + num_states]
        transition_dict[state] = values
    return transition_dict
    
####################################################################################################################
if __name__ == "__main__":
    sample_path = 'ABABBBAAAA'
    states = 'A B'.split()
    extra_dataset_path = 'BBABBBABBAABABABBBAABBBBAAABABABAAAABBBBBAABBABABB'
    print hidden_path_probability(transition_matrix,states,sample_path)