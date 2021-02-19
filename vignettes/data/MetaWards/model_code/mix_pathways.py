from metawards.mixers import merge_matrix_multi_population
from metawards.utils import Console

# use caching to store matrix read in from filename
file_cache = {}

# function to read matrix from file and/or return from cache
def read_file(filename, nage):
    global file_cache
    if filename not in file_cache:
        with open(filename, "r") as FILE:
            contact_matrix = [[num for num in line.split(',')] for line in FILE]
            # set up interaction matrix
            matrix = [[0.0 for i in range(nage)] for j in range(nage)]
            for i in range(nage):
                for j in range(nage):
                    matrix[i][j] = float(contact_matrix[i][j])
        Console.debug("Interaction matrix", variables = [matrix])
        file_cache[filename] = matrix
    return file_cache[filename]
    
# mixing matrix
def mix_pathways(network, **kwargs):
    # extract user parameters
    params = network.params
    
    # extract contact matrix and scaling information
    nage = params.user_params["nage"]
    contact_matrix_filename = params.user_params["contact_matrix_filename"]
    
    # set up interaction matrix in cache or read from cache
    matrix = read_file(contact_matrix_filename, nage)
    
    network.demographics.interaction_matrix = matrix

    return [merge_matrix_multi_population]
