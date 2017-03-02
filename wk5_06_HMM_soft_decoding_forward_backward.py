# -*- coding: utf-8 -*-
'''
Solution to the Soft Decoding Problem.
Finding Mutations in DNA and Proteins (Bioinformatics VI) on Coursera.org
Week 5, code challenge #5
https://stepik.org/lesson/Hidden-Markov-Models-Code-Challenges-(Week-2)-11632/step/12?course=Stepic-Interactive-Text-for-Week-5&unit=9009
'''
import numpy as np

def soft_decoding(transition_matrix,emission_matrix,states,emission_letters,emission_string):
    '''
    Solve the Soft Decoding Problem.
    
    Parameters
    --------
    transition_matrix_input: initial Transition matrix (string)
    emission_matrix_input: initial Emission matrix (string)
    states:  Collection of states of the HMM (list of strings) 
    emission_letters: The HMM's alphabet Σ (list of strings)
    emission_string: A string x of symbols emitted from an HMM (string)
    
    Return
    --------
    None
    
    Print output: An |x| x |States| matrix whose (i, k)-th element holds the conditional probability Pr(πi = k|x).
    '''
    def create_matrices(transition_matrix_input,emission_matrix_input,states,emission_letters):
        '''
        Create numpy arrays from the inputted emission and transition matrices (strings)
        '''
        num_emission_letters = len(emission_letters)
        num_states = len(states)
        transition_matrix = np.zeros((num_states,num_states))
        emission_matrix = np.zeros((num_states,num_emission_letters))
        
        # create the transition matrix
        for i in range(1,num_states+1):
            row = transition_matrix_input.split('\n')[i]
            row = [float(x) for x in row.split()[1:]]
            transition_matrix[i-1] = row

        # create the emission matrix
        for i in range(1,num_states+1):
            row = emission_matrix_input.split('\n')[i]
            row = [float(x) for x in row.split()[1:]]
            emission_matrix[i-1] = row
        return transition_matrix,emission_matrix
        
    transition_matrix, emission_matrix = create_matrices(transition_matrix,emission_matrix,states,emission_letters) 
    reverse_emission_string = emission_string[::-1]
    forward = np.zeros((len(emission_string),len(states)))
    backward = np.zeros((len(emission_string),len(states)))
    
    #create the forward and backward matrices    
    for i in xrange(len(emission_string)):
        if i == 0:
            for j in range(len(states)):
                emission_letter_index = emission_letters.index(emission_string[0])
                forward[0][j] = 1.0/len(states)*emission_matrix[j][emission_letter_index] 
                backward[0][j] = 1.0
            continue
        for j in range(len(states)):
            for k in range(len(states)):
                emission_letter_index = emission_letters.index(emission_string[i])
                reverse_emission_letter_index = emission_letters.index(reverse_emission_string[i-1])
                forward[i][j] += forward[i-1][k]*transition_matrix[k][j]*emission_matrix[j][emission_letter_index]

                backward[i][j] += backward[i-1][k]*transition_matrix[j][k]*emission_matrix[k][reverse_emission_letter_index]
                    
    backward = backward[::-1]
    result_array = np.zeros((len(emission_string),len(states)))

    # combine the forward and backward matrices to create the result matrix
    res = {}
    result_array = np.zeros((len(emission_string),len(states)))
    for i in xrange(len(emission_string)):
        res[i] = {}
        total = 0
        for j in range(len(states)):
            result_array[i][j] = forward[i][j]*backward[i][j]
            total +=  result_array[i][j]
        for j in range(len(states)):
            result_array[i][j] = result_array[i][j]/total
    
    # print the result array
    states_string = ''
    for state in states:
        states_string += state + '\t'
    print states_string
    for i in range(len(emission_string)):
        result_string = ''
        for j in range(len(states)):
            result_string += str(round(result_array[i][j],4)) +'\t'
        print result_string
####################################################################################

if __name__ == "__main__":
    #sample_data
    emission_string = 'zyxxxxyxzz'
    emission_letters = 'x y z'.split()
    states = 'A B'.split()
    
    transition_matrix = '''A	B
    A	0.911	0.089
    B	0.228	0.772'''
    
    emission_matrix = '''x	y	z
    A	0.356	0.191	0.453 
    B	0.040	0.467	0.493'''
    
    soft_decoding(transition_matrix,emission_matrix,states,emission_letters,emission_string)
