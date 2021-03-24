from metawards.mixers import merge_matrix_multi_population
from metawards.utils import Console
from datetime import datetime

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

# determine the lock-down status based on the population and current network
def get_lock_down_status(network, population):

    date = population.date
    params = network.params.user_params

    y1 = int(params["lockdown_date_1_year"])
    m1 = int(params["lockdown_date_1_month"])
    d1 = int(params["lockdown_date_1_day"])

    y2 = int(params["lockdown_date_2_year"])
    m2 = int(params["lockdown_date_2_month"])
    d2 = int(params["lockdown_date_2_day"])

    # Lock down dates
    lock_1 = datetime(y1, m1, d1).date()
    lock_2 = datetime(y2, m2, d2).date()
    
    state = 0
    if date >= lock_1:
        state += 1
    if date >= lock_2:
        state += 1
    return state
    
# mixing matrix
def mix_pathways(network, population, **kwargs):
    # extract user parameters
    params = network.params
    
    # get lockdown status
    state = get_lock_down_status(network = network, population = population)
    
    # extract contact matrix and scaling information
    nage = params.user_params["nage"]
    if state == 0:
        contact_matrix_filename = params.user_params["contact_matrix1_filename"]
    else:
        contact_matrix_filename = params.user_params["contact_matrix2_filename"]
    
    # set up interaction matrix in cache or read from cache
    matrix = read_file(contact_matrix_filename, nage)
    
    network.demographics.interaction_matrix = matrix

    return [merge_matrix_multi_population]
