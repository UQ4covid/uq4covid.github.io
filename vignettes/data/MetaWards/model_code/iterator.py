from metawards.utils import Console
from metawards.iterators import advance_infprob
from metawards.iterators import advance_play
from metawards.iterators import advance_fixed
from metawards.iterators import advance_foi
from metawards.iterators import advance_recovery
from metawards.iterators import advance_foi_work_to_play
from metawards.iterators import advance_work_to_play
from metawards.movers import MoveGenerator, go_ward
from datetime import datetime
import numpy as np
import sys
import re

# use caching to store matrix read in from filename
seed_file_cache = {}

# functions to read seeding probabilities from file and/or return from cache
def read_seed_file(filename):
    global seed_file_cache
    if filename not in seed_file_cache:
        with open(filename, "r") as FILE:
            ward_probs = [[num for num in line.split(',')] for line in FILE]
        Console.debug("Ward-level seeding probabilities:", variables = [ward_probs])
        
        # convert to correct format
        ward_probs = np.array(ward_probs)
        ward_probs_ind = ward_probs[:, 0].astype(int)
        ward_probs_LAD = ward_probs[:, 1].astype(int)
        ward_probs = ward_probs[:, 2].astype(float)
        
        # save in cache        
        seed_file_cache[filename] = "STORED"
        seed_file_cache["ward_probs_ind"] = ward_probs_ind
        seed_file_cache["ward_probs_LAD"] = ward_probs_LAD
        seed_file_cache["ward_probs"] = ward_probs
        
    return seed_file_cache["ward_probs_ind"], seed_file_cache["ward_probs_LAD"], seed_file_cache["ward_probs"]

# read in age lookup
def read_age_file(filename):
    global seed_file_cache
    if filename not in seed_file_cache:
        with open(filename, "r") as FILE:
            age_probs = [[num for num in line.split(',')] for line in FILE]
        Console.debug("Age-level seeding probabilities:", variables = [age_probs])
        
        # convert to correct format
        age_probs = np.array(age_probs)
        age_probs_ind = age_probs[:, 0].astype(int)
        age_probs = age_probs[:, 1].astype(float)
        
        # save in cache
        seed_file_cache[filename] = "STORED"
        seed_file_cache["age_probs_ind"] = age_probs_ind
        seed_file_cache["age_probs"] = age_probs
        
    return seed_file_cache["age_probs_ind"], seed_file_cache["age_probs"]

# functions to read seeding probabilities from file and/or return from cache
def read_states_file(filename):
    global seed_file_cache
    if filename not in seed_file_cache:
        with open(filename, "r") as FILE:
            ini_states = [[num for num in line.split(',')] for line in FILE]
        Console.debug("Initial states:", variables = [ini_states])
        
        # convert to correct format
        ini_states = np.array(ini_states)
        ini_states_output = ini_states[:, 0]
        ini_states_LAD = ini_states[:, 1].astype(int)
        ini_states = ini_states[:, 2:].astype(int)
        
        # save in cache        
        seed_file_cache[filename] = "STORED"
        seed_file_cache["ini_states_output"] = ini_states_output
        seed_file_cache["ini_states_LAD"] = ini_states_LAD
        seed_file_cache["ini_states"] = ini_states
        
    return seed_file_cache["ini_states_output"], seed_file_cache["ini_states_LAD"], seed_file_cache["ini_states"]
    
# determine the lock-down status based on the population and current network
def get_lock_down_vars(network, population):

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
    rate = 1.0
    if date >= lock_1:
        state += 1
        rate = params["lock_1_restrict"]
    if date >= lock_2:
        state += 1
        rate = (1.0 - ((1.0 - params["lock_1_restrict"]) * params["lock_2_release"]))

    can_work = params["can_work"][state]
    return state, rate, can_work

# advance in a lock-down state weekday
def advance_lockdown_week(network, population, **kwargs):

    state, rate, can_work = get_lock_down_vars(network = network, population = population)

    advance_infprob(scale_rate = rate, network = network, population = population, **kwargs)
    advance_play(network = network, population = population, **kwargs)
    if can_work:
        advance_fixed(network = network, population = population, **kwargs)

# advance in a lockdown state weekend        
def advance_lockdown_weekend(network, population, **kwargs):

    state, rate, can_work = get_lock_down_vars(network = network, population = population)
    
    ## update dyn_play_at_home
    dyn_play_at_home = []
    for i in range(0, len(network.subnets)):
        dyn_play_at_home.append(network.subnets[i].params.dyn_play_at_home)
        network.subnets[i].params.dyn_play_at_home = network.subnets[i].params.user_params["p_home_weekend"]

    advance_infprob(scale_rate = rate, network = network, population = population, **kwargs)
    advance_work_to_play(network = network, population = population, **kwargs)
    
    ## reset dyn_play_at_home
    for i in range(0, len(network.subnets)):
        network.subnets[i].params.dyn_play_at_home = dyn_play_at_home[i]
        
# advance FOI weekend
def advance_foi_work_to_play_weekend(network, population, **kwargs):
    
    ## update dyn_play_at_home
    dyn_play_at_home = []
    for i in range(0, len(network.subnets)):
        dyn_play_at_home.append(network.subnets[i].params.dyn_play_at_home)
        network.subnets[i].params.dyn_play_at_home = network.subnets[i].params.user_params["p_home_weekend"]

    advance_foi_work_to_play(network = network, population = population, **kwargs)
    
    ## reset dyn_play_at_home
    for i in range(0, len(network.subnets)):
        network.subnets[i].params.dyn_play_at_home = dyn_play_at_home[i]

# set custom advance function
def advance_initial_seeds(output_dir, network, population, infections, profiler, rngs, **kwargs):
    
    # extract user parameters
    params = network.params
    
    # extract files name for initial seeding probabilities
    ward_seed_filename = params.user_params["ward_seed_filename"]
    age_seed_filename = params.user_params["age_seed_filename"]
    ini_states_filename = params.user_params["ini_states_filename"]
    
    # start profiler
    p = profiler.start("additional_seeds")
    
    # set up lookups or read from cache
    age_probs_ind, age_probs = read_age_file(age_seed_filename)
    ward_probs_ind, ward_probs_LAD, ward_probs = read_seed_file(ward_seed_filename)
    ini_states_output, ini_states_LAD, ini_states = read_states_file(ini_states_filename)
    
    # get output identifier for this run of MetaWards
    output = re.search(r"Ens([0-9a-z]*)", output_dir.output_dir())
    output = output.group(0)
    
    # set up vector of state names
    state_names = ["E", "P", "I1", "I2", "RI", "DI", "A", "RA", "H", "RH", "DH"]
    state_names = np.array(state_names)
    state_nums = range(len(state_names))
    state_nums = np.array(state_nums)
    
    # extract states for this run
    filter_states = [i == output for i in ini_states_output]
    if(sum(filter_states) == 0):
        raise Exception(f"Cannot find {output} seeding information.")
    
    # filter seeding information relating to this run
    ini_states_LAD = ini_states_LAD[filter_states]
    ini_states = ini_states[filter_states, :]
    
    # loop over LADs
    for k in range(len(ini_states_LAD)):
        
        # check for any seeds
        filter_LAD = [i == ini_states_LAD[k] for i in ini_states_LAD]
        ini_states_curr = ini_states[filter_LAD, :][0]
    
        # loop over states
        for j in range(len(ini_states_curr)):
    
            # extract number of affected individuals
            nind = ini_states_curr.sum()
            
            if nind > 0:
                
                # extract wards
                filter_wards = [i == ini_states_LAD[k] for i in ward_probs_LAD]
                if sum(filter_wards) == 0:
                    raise Exception(f"Can't find any wards")
                tward_probs = ward_probs[filter_wards]
                tward_probs_ind = ward_probs_ind[filter_wards]
                
                # select seeds in age-classes at random according to initial probabilities
                age_seeds = np.random.multinomial(nind, age_probs)
                
                # run over each age demographic
                for demographic in range(len(age_seeds)):
                    
                    # check if any seeding done in demographic
                    while age_seeds[demographic] > 0:
                        
                        # select states for given seeds
                        state_probs = ini_states_curr / ini_states_curr.sum()
                        
                        # select states
                        states = np.random.multinomial(1, state_probs)
                        filter_state = [i == 1 for i in states]
                        state = state_names[filter_state]
                        staten = state_nums[filter_state]
                        
                        # remove selected state from states
                        ini_states_curr[staten] -= 1
                        
                        # select seeds in wards at random according to initial probabilities
                        ward = np.random.multinomial(1, tward_probs)
                        filter_ward = [i == 1 for i in ward]
                        ward = tward_probs_ind[filter_ward]
                        
                        # generate move
                        move = MoveGenerator(from_demographic = f'age{demographic + 1}',
                            to_demographic = f'age{demographic + 1}',
                            from_ward = ward,
                            to_ward = ward,
                            from_stage = "S",
                            to_stage = state[0],
                            number = 1)

                        go_ward(generator = move, output_dir = output_dir, network = network, 
                                population = population, infections = infections, 
                                profiler = profiler, rngs = rngs, **kwargs)
                            
                        # update counter
                        age_seeds[demographic] -= 1

    # end profiler
    p.stop()

# custom iterator function
def iterate_lockdown(stage: str, population, **kwargs):

    # is this a weekday or a weekend?
    if population.date is None:
        # have to guess from the day - assume day zero was a Monday
        day = population.day % 7
        # 0-4 is a weekday, 5 and 6 are weekend
        is_weekend = (day >= 5)
    else:
        day = population.date.weekday()
        # 0-4 is a weekday, 5 and 6 are weekend
        is_weekend = (day >= 5)

    # need custom advance functions to scale rates after lockdown
    # in period before lockdown this treats workers as players at
    # weekends
    if population.day == 1:
        if stage == "setup":
            return [advance_initial_seeds]
        else:
            return []
    else:
        if stage == "foi":
            if is_weekend:
                return [advance_foi_work_to_play_weekend, advance_recovery]
            else:
                return [advance_foi, advance_recovery]
        elif stage == "infect":
            if is_weekend:
                return [advance_lockdown_weekend]
            else:
                return [advance_lockdown_week]
        else:
            # we don't do anything at the "analyse" or "finalise" stages
            return []

