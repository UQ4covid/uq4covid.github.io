from metawards.mixers import merge_using_matrix
from metawards.utils import Console

# use caching to store matrix read in from filename
file_cache = {}

# function to read matrix from file and/or return from cache
def read_file(filename, nage, nu, betaA_GP):
    global file_cache
    if filename not in file_cache:
        with open(filename, "r") as FILE:
            contact_matrix = [[num for num in line.split(',')] for line in FILE]
            # set up interaction matrix
            matrix = [[0.0 for i in range(nage * 3)] for j in range(nage * 3)]
            for i in range(nage):
                for j in range(nage):
                    matrix[i][j] = float(contact_matrix[i][j]) * nu
            for i in range(nage):
                for j in range(nage):
                    matrix[i][j + nage] = matrix[i][j] * betaA_GP
        Console.debug("Interaction matrix", variables = [matrix])
        file_cache[filename] = matrix
    return file_cache[filename]
    
# mixing matrix
def mix_pathways(network, **kwargs):
    # extract user parameters
    params = network.params
    
    # extract contact matrix and scaling information
    nage = params.user_params["nage"]
    nu = params.user_params["nu"]
    contact_matrix_filename = params.user_params["contact_matrix_filename"]
    
    # extract how much FOI to GP population is affected by others
    betaA_GP = params.user_params["betaA_GP"]
    
    # set up interaction matrix in cache or read from cache
    matrix = read_file(contact_matrix_filename, nage, nu, betaA_GP)
    
    network.demographics.interaction_matrix = matrix

    return [merge_using_matrix]
