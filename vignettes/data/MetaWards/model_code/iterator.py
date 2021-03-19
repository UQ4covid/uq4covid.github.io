from metawards.utils import Console
from metawards.iterators import advance_infprob
from metawards.iterators import advance_play
from metawards.iterators import advance_fixed
from datetime import datetime
import numpy as np

# use caching to store matrix read in from filename
seed_file_cache = {}

# functions to read seeding probabilities from file and/or return from cache
def read_seed_file(filename):
    global seed_file_cache
    if filename not in seed_file_cache:
        with open(filename, "r") as FILE:
            ward_probs = [[num for num in line.split(',')] for line in FILE]
        Console.debug("Ward-level seeding probabilities:", variables = [ward_probs])
        seed_file_cache[filename] = ward_probs
    return seed_file_cache[filename]
    
def read_age_file(filename):
    global seed_file_cache
    if filename not in seed_file_cache:
        with open(filename, "r") as FILE:
            age_probs = [[num for num in line.split(',')] for line in FILE]
        Console.debug("Age-level seeding probabilities:", variables = [age_probs])
        seed_file_cache[filename] = age_probs
    return seed_file_cache[filename]
    
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

    advance_infprob(scale_rate = rate, network = network, population = population, **kwargs)
    advance_play(network = network, population = population, **kwargs)

# set custom advance function
def advance_initial_seeds(network, population, infections, profiler, rngs, **kwargs):
    
    # extract user parameters
    params = network.params
    
    # extract files name for initial seeding probabilities
    ward_seed_filename = params.user_params["ward_seed_filename"]
    age_seed_filename = params.user_params["age_seed_filename"]
    
    # start profiler
    p = profiler.start("additional_seeds")
    
    # only generate seeds once
    if ward_seed_filename and age_seed_filename not in seed_file_cache:
        # extract required number of seeds
        nseeds = params.user_params["nseeds"]
        
        # set up age-class probabilities in cache or read from cache
        age_probs = read_age_file(age_seed_filename)
        age_probs = np.array(age_probs)
        age_probs_ind = age_probs[:, 0].astype(int)
        age_probs = age_probs[:, 1].astype(float)
        
        # select seeds in age-classes at random according to initial probabilities
        age_seeds = np.random.multinomial(nseeds, age_probs)
        
        # set up ward probabilities in cache or read from cache
        ward_probs = read_seed_file(ward_seed_filename)
        ward_probs = np.array(ward_probs)
        ward_probs_ind = ward_probs[:, 0].astype(int)
        ward_probs = ward_probs[:, 1].astype(float)
        
        # run over each age demographic
        for demographic in range(len(age_seeds)):
            
            ## check if any seeding done in demographic
            if age_seeds[demographic] > 0:
                
                # select seeds in wards at random according to initial probabilities
                seeds = np.random.multinomial(age_seeds[demographic], ward_probs)
                    
                # now seed infections
                for i in range(len(seeds)):
                    ward = ward_probs_ind[i]
                    num = seeds[i]
                    
                    if num > 0:
                        seed_network = network.subnets[demographic]
                        seed_wards = seed_network.nodes
                        seed_infections = infections.subinfs[demographic].play
                        
                        try:
                            ward = seed_network.get_node_index(ward)
            
                            if seed_wards.play_suscept[ward] == 0:
                                Console.warning(
                                    f"Cannot seed {num} infection(s) in ward {ward} "
                                    f"as there are no susceptibles remaining")
                                continue
            
                            elif seed_wards.play_suscept[ward] < num:
                                Console.warning(
                                    f"Not enough susceptibles in ward to see all {num}")
                                num = seed_wards.play_suscept[ward]
            
                            seed_wards.play_suscept[ward] -= num
                            if demographic is not None:
                                Console.print(
                                    f"seeding demographic {demographic} "
                                    f"play_infections[0][{ward}] += {num}")
                            else:
                                Console.print(
                                    f"seeding play_infections[0][{ward}] += {num}")
            
                            seed_infections[0][ward] += num
            
                        except Exception as e:
                            Console.error(
                                f"Unable to seed the infection using {seed}. The "
                                f"error was {e.__class__}: {e}. Please double-check "
                                f"that you are trying to seed a node that exists "
                                f"in this network.")
                            raise e

    # end profiler
    p.stop()

# custom iterator function
def iterate_lockdown(population, **kwargs):

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
    # in period before lockdown this is the same as iterate_working_week
    if is_weekend:
        return [advance_initial_seeds, advance_lockdown_weekend]
    else:
        return [advance_initial_seeds, advance_lockdown_week]
