"""
Solution to the Burrows-Wheeler Transform Construction Problem.
Finding Mutations in DNA and Proteins (Bioinformatics VI) on Coursera.org
Week 2, code challenge #2
https://stepik.org/lesson/The-Burrows-Wheeler-Transform-297/step/5?unit=8999
"""

def burrows_wheeler(seq):
    '''
    Construct the Burrows-Wheeler transform of a given seq.
    
    Parameters
    --------
    seq: A DNA sequence
    
    Return
    --------
    Burrow-wheeler transform of the input sequence (last column of the array) (string)
    '''
    bw_array = [seq]# burrows_wheeler array
    bw_seq = ''
    for i in range(1,len(seq)):
        shifted_seq = seq[i:] + seq[:i]
        bw_array.append(shifted_seq)
    bw_array = sorted(bw_array)
    for shifted_seq in bw_array:
        bw_seq += shifted_seq[-1]
    return bw_seq
    
    
#################################################
if __name__ == "__main__":
    sample_input = 'GCGTGCCTGGTCA$'
    bw_seq = burrows_wheeler(sample_input)
    print bw_seq


