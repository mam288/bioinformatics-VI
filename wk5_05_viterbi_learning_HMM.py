'''
Solution to Implement Viterbi Learning
Finding Mutations in DNA and Proteins (Bioinformatics VI) on Coursera.org
Week 5,code challenge #5
https://stepik.org/lesson/Hidden-Markov-Models-Code-Challenges-(Week-2)-11632/step/10?course=Stepic-Interactive-Text-for-Week-5&unit=9009 
'''
import networkx as nx
import numpy as np

def viterbi_learn(input_transition_matrix,input_emission_matrix,states,emission_letters,emission_string,iterations):
    '''
    Implement Viterbi learning for estimating the parameters of an HMM.
    
    Parameters
    --------
    input_transition_matrix: Initial transition matrix for the HMM (string)
    input_emission_matrix: Initial emission matrix for the HMM (string)
    states: the HMM's states (list of strings)
    emission_letters: HMM's alphabet Σ (list of strings)
    emission_string: string x of symbols emitted by an HMM (string)
    iterations: a number of iterations j (integer)

    Returns
    --------
    None 
    
    Print output: Emission and transition matrices resulting from applying 
    Viterbi learning for j iterations.
    '''
    hidden_path = viterbi(input_transition_matrix,input_emission_matrix,states,emission_letters,emission_string) # find the hidden path for the givven matrices and emission string
    new_emission_matrix,new_transition_matrix = estimate_HMM_parameters(emission_string,hidden_path,states,emission_letters) # use the path to calculate new emission and transition matrices
    
    # for a given number of iterations use the new matrices to calculate a new hidden path and then use that path to calculate new matrices
    for i in range(iterations):
        hidden_path = viterbi(new_transition_matrix,new_emission_matrix,states,emission_letters,emission_string)
        new_emission_matrix,new_transition_matrix = estimate_HMM_parameters(emission_string,hidden_path,states,emission_letters)
    
    print_matrix(new_transition_matrix,states,states)
    print ('--------')
    print_matrix(new_emission_matrix,states,emission_letters)
    
def estimate_HMM_parameters(emission_string,hidden_path,states,emission_letters):
    '''
    Estimate the transition matrix and emission matrix given a hidden path and emission string
    '''
    def create_transition_matrix(emission_string,hidden_path,states):
        '''
        Creates a transition matrix using the hidden path,emission string and states.
        Output: Transition matrix as a numpy array 
        '''
        transition_matrix = np.zeros((len(states),len(states)))
        
        # add up the number of times the hidden path transitions from each state to each other state
        for i in range(len(hidden_path)-1):
            state_1 = hidden_path[i]; state_2 = hidden_path[i+1]   # get the state at i and i+1 in the hidden path
            state_1_index = states.index(state_1); state_2_index = states.index(state_2)  # get the index for each state
            transition_matrix[state_1_index][state_2_index] += 1

        # add pseudo counts and normalize the matrix
        for i in range(len(states)):
            row_sum = sum(transition_matrix[i])
            if row_sum == 0:
                transition_matrix[i] += 0.01  # add pseudo counts
            transition_matrix[i] = transition_matrix[i]/row_sum
        return transition_matrix
        
    def create_emission_matrix(emission_string,hidden_path,states,emission_letters):
        emission_matrix = np.zeros((len(states),len(emission_letters)))
        
        # add up the number of times a letter is emitted for each state in the hidden path
        for i in range(len(hidden_path)):
            emission = emission_string[i]
            emission_index = emission_letters.index(emission)
            state = hidden_path[i]
            state_index = states.index(state)
            emission_matrix[state_index][emission_index] +=1    
    
        # add pseudo counts and normalize the matrix
        for i in range(len(states)):
            row_sum = sum(emission_matrix[i])
            if row_sum == 0:
                emission_matrix[i] += 0.01
            emission_matrix[i] = emission_matrix[i]/row_sum
        return emission_matrix
    emission_matrix = create_emission_matrix(emission_string,hidden_path,states,emission_letters)
    transition_matrix = create_transition_matrix(emission_string,hidden_path,states)
    return emission_matrix,transition_matrix

def viterbi(transition_matrix_input,emission_matrix_input,states,emission_letters,emission_string):
    def create_matrices(transition_matrix_input,emission_matrix_input,states,emission_letters):
        '''
        Convert the transition and emission matrices from strings to numpy arrays
        '''
        num_emission_letters = len(emission_letters)
        num_states = len(states)
        transition_matrix = np.zeros((num_states,num_states))
        emission_matrix = np.zeros((num_states,num_emission_letters))
        for i in range(1,num_states+1):
            t_row = transition_matrix_input.split('\n')[i]
            t_row = [float(x) for x in t_row.split()[1:]]
            transition_matrix[i-1] = t_row
            e_row = emission_matrix_input.split('\n')[i]
            e_row = [float(x) for x in e_row.split()[1:]]
            emission_matrix[i-1] = e_row
        return emission_matrix,transition_matrix
    
    def create_viterbi_graph(transition_matrix_input,emission_matrix_input,states,emission_letters,emission_string):
        '''
        Create the viterbi graph with weighted edges using the given transition and emission matrices and emission string.
        '''
        viterbi_graph = nx.DiGraph()
        num_states = len(states)
        if type(transition_matrix_input) == np.ndarray:
            transition_matrix = transition_matrix_input
            emission_matrix = emission_matrix_input
        else:
            # convert the matrices to numpy arrays if they are not already arrays.
            emission_matrix,transition_matrix = create_matrices(transition_matrix_input,emission_matrix_input,states,emission_letters)
        
        for i in range(len(emission_string)+1):
            if i == len(emission_string):
                index = 0
            else:  # find the index of the letter emitted at i
                index = emission_letters.index(emission_string[i])
                
            # calculate the 'weight' for each edge in the graph
            for outgoing_index in range(len(states)):
                outgoing_state =   states[outgoing_index]
                emission_weight = emission_matrix[outgoing_index][index]
                if i == 0:
                    incoming_state = ''
                    transition_weight = 0.5
                    viterbi_graph.add_edge(str(i),str(i+1) + outgoing_state,{'weight': emission_weight*transition_weight})
                else:
                    for incoming_index in range(num_states):
                        incoming_state = states[incoming_index]
                        if i == len(emission_string):
                            transition_weight = 1; emission_weight = 1; outgoing_state = ''
                            viterbi_graph.add_edge(str(i)  + incoming_state,str(i+1),{'weight': emission_weight*transition_weight})
                            continue 
                        transition_weight = transition_matrix[incoming_index][outgoing_index]                     
                        viterbi_graph.add_edge(str(i) + incoming_state,str(i+1) + outgoing_state,{'weight': emission_weight*transition_weight})
        return viterbi_graph
    
    def calculate_scores_viterbi(viterbi_graph):
        '''
        Go through the viterbi graph in topological order and score each node.
        '''
        viterbi_graph.node['0']['score'] = 1
        sorted_nodes = nx.topological_sort(viterbi_graph)
        backtrack = {}
        for node in sorted_nodes[1:]:
            predecessors = viterbi_graph.predecessors(node)
            high_score = 0
            backtrack_node = ''
            
            # calculate the score from each predecessor to 'node' and set it as the high score if it is greater than the current high score
            for predecessor in predecessors:
                score = viterbi_graph.node[predecessor]['score']*viterbi_graph[predecessor][node]['weight']
                if score > high_score:
                    high_score = score
                    backtrack_node = predecessor
                    
            viterbi_graph.node[node]['score'] = high_score
            backtrack[node] = backtrack_node # record which predecessor was used to calculate the high score
        return backtrack
        
    def reconstruct_path(viterbi_graph,backtrack):
        '''Reconstruct the hidden path from the scored viterbi_graph'''
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
    backtrack= calculate_scores_viterbi(viterbi_graph)
    path = reconstruct_path(viterbi_graph,backtrack)
    return path
    
def print_matrix(matrix,rows,columns):
    '''Print the inputted matrix'''
    column_string = '\t'
    for r in range(len(columns)):
        column_string += columns[r] + '\t'
    print column_string[:-1]
    for r in range(len(rows)):
        row = rows[r] + '\t'
        for c in range(len(columns)):
            num = matrix[r][c]
            num = round(num,3)
            row += str(num) + '\t'
        print row[:-1]
    
####################################################################################
if __name__ == "__main__":
    #sample_data
    transition_matrix = '''A	B
    A	0.599	0.401	
    B	0.294	0.706'''	
    emission_matrix = '''x	y	z
    A	0.424	0.367	0.209	
    B	0.262	0.449	0.289'''
    iterations = 10
    emission_string = 'zyzxzxxxzz'
    emission_letters = 'x y z'.split()
    states = 'A B'.split()
    
    viterbi_learn(transition_matrix,emission_matrix,states,emission_letters,emission_string,iterations)