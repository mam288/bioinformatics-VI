# -*- coding: utf-8 -*-
'''
Solution to HMM Parameter Estimation Problem.
Finding Mutations in DNA and Proteins (Bioinformatics VI) on Coursera.org
Week 5, code challenge #4 
https://stepik.org/lesson/Hidden-Markov-Models-Code-Challenges-(Week-2)-11632/step/8?course=Stepic-Interactive-Text-for-Week-5&unit=9009
'''
import numpy as np

def estimate_HMM_parameters(emission_string,hidden_path,states,emission_letters):
    '''
    Solve the HMM Parameter Estimation Problem.
    
    Parameters
    ----------
    emission_string: A string x of symbols emitted from an HMM (string)
    hidden_path: path π (list of strings)
    states:  Collection of states of the HMM (list of strings)      
    emission_letters: The HMM's alphabet Σ (list of strings)
        
    Returns
    -------
    None

    Print output: A transition matrix Transition followed by an emission matrix
        Emission that maximize Pr(x, π) over all possible transition and
        emission matrices.
    '''
    def create_transition_matrix(emission_string, hidden_path, states):
        '''
        Create a transition matrix using the given emission string, hidden path,
        and states. 
        '''
        transition_matrix = np.zeros((len(states), len(states)))

        # add up the number of times the hidden path transitions from each state to each other state
        for i in range(len(hidden_path)-1):
            state_1 = hidden_path[i]; state_2 = hidden_path[i+1]
            state_1_index = states.index(state_1); state_2_index = states.index(state_2)
            transition_matrix[state_1_index][state_2_index] += 1

        # add pseudo counts and normaize the matrix
        for i in range(len(states)):
            row_sum = sum(transition_matrix[i])
            if row_sum == 0:
                transition_matrix[i] += 0.01
                row_sum = sum(transition_matrix[i])
            transition_matrix[i] = transition_matrix[i]/row_sum
        print_matrix(transition_matrix,states,states)

    def create_emission_matrix(emission_string,hidden_path,states,emission_letters):
        '''
        Create an emission  matrix using the given emission string, hidden path,
        emission_leters and states.
        '''
        emission_matrix = np.zeros((len(states),len(emission_letters)))

        # add up the number of times a letter is emitted for each state in the hidden path
        for i in range(len(hidden_path)):
            emission = emission_string[i] # find out which letter was emitted at i
            emission_index = emission_letters.index(emission)
            state = hidden_path[i]   # find out the state at i
            state_index = states.index(state)
            emission_matrix[state_index][emission_index] +=1

        # add pseudo counts normalize the matrix
        for i in range(len(states)):
            row_sum = sum(emission_matrix[i])
            if row_sum == 0:
                emission_matrix[i] += 0.01
                row_sum = sum(emission_matrix[i])
            emission_matrix[i] = emission_matrix[i]/row_sum

        print_matrix(emission_matrix,states,emission_letters)

    create_transition_matrix(emission_string,hidden_path,states)
    print ('--------')
    create_emission_matrix(emission_string,hidden_path,states,emission_letters)

def print_matrix(matrix,rows,columns):
    '''
    Prints a given emission or transition matrix.
    '''
    column_string = '\t'
    for r in range(len(columns)):
        column_string += columns[r] + '\t'
    print column_string[:-1]
    for r in range(len(rows)):
        row = rows[r] + '\t'
        for c in range(len(columns)):
            num = matrix[r][c]
            row += str(round(num,3)) + '\t'  # round the number to the 3rd decimal place
        print row[:-1]

###################################################################################
if __name__ == "__main__":
    #sample_data
    emission_string = 'yzzzyxzxxx'
    emission_letters = 'x y z'.split()
    hidden_path = 'BBABABABAB'
    states = 'A B C'.split()

    estimate_HMM_parameters(emission_string,hidden_path,states,emission_letters)