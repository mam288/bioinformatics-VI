'''
Solution to the Profile HMM with Pseudocounts Problem.
Finding Mutations in DNA and Proteins (Bioinformatics VI) on Coursera.org
Week 5, code challenge #2
https://stepik.org/lesson/Hidden-Markov-Models-Code-Challenges-(Week-2)-11632/step/4?course=Stepic-Interactive-Text-for-Week-5&unit=9009
'''


import networkx as nx
import numpy as np

   
def create_HMM_graph_pseudocounts(emission_letters,alignment,t,ps):
    '''
    Profile HMM with Pseudocounts Problem.
    
    Parameters
    --------
    emission_letters: alphabet Σ
    alignment: multiple alignment Alignment whose strings are formed from Σ
    t: threshold θ
    
    Return
    -------
    None
        
    Print output: The transition and emission matrices of HMM(Alignment, θ, σ).
    
    Note: A viterbi graph with weighted edges is created along with the matrices. 
    '''
    def create_profile_counts():  # t= threshold
        '''
        Count the number of each emmitted letter for each column in the alignment.
        '''
        emission_counts = {}
        for index in range(align_len):
            emission_counts[index+1] = {}
            for letter in emission_letters_w_dash:
                emission_counts[index +1][letter] = 0
            for seq_num in range(num_align):
                emission = alignment[seq_num][index]
                emission_counts[index+1][emission] += 1
        return emission_counts
        
    def find_next_column(seq,index,column_map):
        '''
        Find the column corresponding to the next location in seq that is populated (not a '-'). Returns the column index 
        and a 0 or 1 depending on whether the column is the first or second entry at that location in the column map. 
        '''
        for i in range(index,len(seq)):
            column = column_map[i+1]
            char = seq[i]
            if char != '-': 
                return i,0
            elif char == '-' and len(column) >1:
                return i,1
        return len(seq),0
        
        
    def transition_graph(graph):
        '''Create a transition matrix and add weights to the edges in the viterbi graph.'''
        transition_matrix = np.zeros((len(node_list),len(node_list)))
        matrix_column_map = {}
        for i,node in enumerate(node_list):
            matrix_column_map[node] = i
        matrix_column_map 
        adjustment_factor = 1
        
        # for each node, add outgoing edges and successor nodes to the graph
        for node in graph.nodes():
            successor_nodes=[]
            if len(node) >1:
                index = int(node[1:])
            else:
                index = 0
            if node_list.index(node) > len(node_list) -5: #if the node is in the last 4 columns
                successor_nodes += ['E']   # add 'E' (sink node) to the list of successors
                
            # create names for successor nodes
            if node[0] == 'I':
                successor_nodes += [node,'M'+str(index+1),'D'+str(index+1)] 
            else:
                successor_nodes += ['M'+str(index+1),'D'+str(index+1),'I' + str(index)]
                                    
            # initialize successor nodes and transition matrix with the pseudocount value
            for successor in successor_nodes:
                if successor in graph.nodes() and node != 'E':
                    graph.add_edge(node,successor,weight = ps)
                    transition_matrix[node_list.index(node)][node_list.index(successor)] = ps
                                      
        for n in range(align_len+1):
            char_index = n-1
            current_column_list = column_map[n]
            adjustment_factor = column_map.values().count(current_column_list)
            for seq_index in range(num_align):
                current_char = alignment[seq_index][char_index]

                # if the column is a 'M/D' column
                if len(current_column_list) >1:
                    if current_char != '-':
                        current_column  = current_column_list[0] #  choose 'M' column
                    elif current_char == '-':
                        current_column = current_column_list[1]    # choose 'D' column

                # if the column is an 'I' column           
                else:
                    current_column = current_column_list[0]
                    if current_char == '-' and current_column[0] == 'I':
                        continue
              
                next_column_index,list_index = find_next_column(alignment[seq_index],n,column_map)
                next_column = column_map[next_column_index+1][list_index]

                # calculate the number of dashes is the column
                if current_column == 'S':
                    dash_count = 0
                    char_count = num_align
                else:
                    dash_count = emission_counts[n]['-']
                    char_count = (num_align - dash_count)
                dash_count = dash_count*adjustment_factor
                char_count = char_count*adjustment_factor

                # calculate the transition probabilities and populate the transition matrix and graph
                if current_column[0] == 'D':
                    transition_matrix[matrix_column_map[current_column]][matrix_column_map[next_column]] += 1/float(dash_count) 
                    graph[current_column][next_column]['weight'] += 1/float(dash_count) 
                else:
                    transition_matrix[matrix_column_map[current_column]][matrix_column_map[next_column]] += 1/float(char_count)
                    graph[current_column][next_column]['weight'] += 1/float(char_count) 
                    
        # normalize the transition matrix and graph probabilities
        for i in range(len(transition_matrix[:-1])):
            row_sum = sum(transition_matrix[i])
            successor_list =  graph.successors(node_list[i])
            if row_sum != 0:
                transition_matrix[i] = transition_matrix[i]/row_sum
            for successor in successor_list:
                if row_sum != 0:
                    graph.edge[node_list[i]][successor]['weight']/= row_sum
        return transition_matrix

            
    def create_emission_matrix():
        '''Create an emission matrix and create a viterbi graph.'''
        def reverse_map(column_map):
            '''Reverse the keys and values for the column map'''
            reversed_column_map = {}
            for column_index in column_map:
                for column_name in column_map[column_index]:
                    reversed_column_map[column_name] = column_index
            return reversed_column_map
            
        graph = nx.DiGraph()
        node_list = ['S','I0']
        graph.add_nodes_from(node_list,emission_probabilities =0)
        emission_matrix = np.zeros((10,len(emission_letters)))
        emission_matrix[node_list.index('I0')] += ps
        row_num = 1
        column_map = {0:['S']}


        for i in range(1,align_len+1):
            num_dashes = emission_counts[i]['-']
            dash_percent = num_dashes/float(num_align)
            
            # populate the emission matrix with counts and add nodes to the graph
            if dash_percent >= t:   # if the number of dashes exceeds the given threshold
                column_map[i] = ['I' + str(row_num-1)]
                insertion_row = 'I' + str(row_num-1) 
                graph.add_nodes_from([insertion_row],emission_probabilities =0)
                if insertion_row not in node_list:
                    emission_matrix[node_list.index(insertion_row)] = np.zeros(len(emission_letters))
                    emission_matrix[node_list.index(insertion_row)] += ps
                for letter in range(len(emission_letters)):
                    emission_matrix[node_list.index(insertion_row)][letter] += emission_counts[i][emission_letters[letter]]
            else:
                column_map[i] = ['M' + str(row_num),'D' + str(row_num)]
                rows = ['M' + str(row_num),'D' + str(row_num),'I' + str(row_num)]
                node_list += rows
                graph.add_nodes_from(node_list,emission_probabilities = np.zeros(len(emission_letters)))
                for row in rows:
                    for letter in range(len(emission_letters)):
                        emission_matrix[node_list.index(rows[0])][letter] = emission_counts[i][emission_letters[letter]]
                row_num += 1
 
        # calculate the probabilities for the emission matrix and graph
        total_populated_rows_dict = {}
        for i in range(1,align_len+1):
            rows = column_map[i]
            num_dashes = emission_counts[i]['-']
            total_populated_rows = num_align-num_dashes
            for row in rows: 
                if row in total_populated_rows_dict:
                    total_populated_rows_dict[row] += total_populated_rows
                else:
                    total_populated_rows_dict[row] = total_populated_rows
        rows = total_populated_rows_dict.keys()
        
        # normalize the emission matrix and graph probabilities
        for row in node_list:
            if row[0] == 'D' or row[0] == 'S' or row[0] == 'E':
                continue
            num_dashes = emission_counts[i]['-']
            total_populated_rows = num_align-num_dashes          
            if row in total_populated_rows_dict:
                emission_matrix[node_list.index(row)] = emission_matrix[node_list.index(row)]/float(total_populated_rows_dict[row])
            emission_matrix[node_list.index(row)] += ps
            row_sum = sum(emission_matrix[node_list.index(row)])
            emission_matrix[node_list.index(row)] = emission_matrix[node_list.index(row)]/row_sum
            graph.node[row]['emission_probabilities'] = emission_matrix[node_list.index(row)]

        node_list += ['E']
        column_map[align_len+1] = ['E']
        emission_matrix[node_list.index('E')] = np.zeros(len(emission_letters))
        graph.add_nodes_from(['E'],emission_probabilities = np.zeros(len(emission_letters)))
        return emission_matrix,node_list,column_map,graph
        
    alignment = alignment.split()
    emission_letters_w_dash = emission_letters + ['-']
    align_len = len(alignment[0]) # length of each alignment sequence
    num_align = len(alignment)  # number of sequences in the alignment
    graph = nx.DiGraph
    emission_counts = create_profile_counts()
    emission_matrix,node_list,column_map,graph = create_emission_matrix()
    transition_matrix = transition_graph(graph)
    output_data_HMM_profile(transition_matrix,node_list,node_list)
    print ('--------')
    output_data_HMM_profile(emission_matrix,node_list,emission_letters)
   
def output_data_HMM_profile(matrix,rows,columns):
    '''Print the given matrix.'''
    column_string = '\t'
    for r in range(len(columns)):
        column_string += columns[r] + '\t'

    print column_string[:-1]
    for r in range(len(rows)):
        row = rows[r] + '\t'
        for c in range(len(columns)):
            num = matrix[r][c]
            if num == 0:
                num = 0
            else:
                num = round(num,3)
            row += str(num) + '\t'
        print row[:-1]
####################################################################################
if __name__ == "__main__":
    t = 0.358
    emission_letters = 'A B C D E'.split()
    alignment = '''A-A
    ADA
    ACA
    A-C
    -EA
    D-A'''
    ps = 0.01
    
    create_HMM_graph_pseudocounts(emission_letters,alignment,t,ps)