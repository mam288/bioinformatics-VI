"""
Solution to Suffix Array Construction Problem.
Finding Mutations in DNA and Proteins (Bioinformatics VI) on Coursera.org
Week 2, code challenge #1
https://stepik.org/lesson/Suffix-Arrays-310/step/2?course=Stepic-Interactive-Text-for-Week-2&unit=8998
"""

def generate_suffix_array(seq):
    '''
    Construct the suffix array of a string.
    
    Parameters
    --------
    seq: DNA sequence (string)
    
    Return
    --------
    Suffix array - a string of integers with each integer being 
    the starting point of a suffix. The array is sorted alphabetically by suffix. (string)
    '''
    suffixes = sorted(generate_suffixes_tuples(seq))
    index_string = ''
    for suffix in suffixes:
        index_string += str(suffix[1]) + ', '
    return index_string
    
    
def generate_suffixes_tuples(pattern):
    '''Generates all of the suffixes of pattern.'''
    suffixes = []
    index = 0
    while pattern != '':
        suffixes.append((pattern,index))
        index += 1
        pattern = pattern[1:]
    return suffixes
   
#########################################################################################
if __name__ == "__main__":
    sample_input = 'AACGATAGCGGTAGA$'
    print generate_suffix_array(sample_input)