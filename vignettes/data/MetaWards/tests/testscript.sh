#!/bin/bash

## first test checks a simple SEPID model in a single population
## using the standard MW approach specifying parameters
## through the beta[] and progress[] rates, and setting
## no age-specific mixing and seeding only in one demographic
## these outbeaks can then be compared to deterministic and
## stochastic models developed in R

cd test1
./testscript.sh
cd ..

## same as test1 except except with age-specific mixing thus 
## same R0 = 3 results in different nu. Only discrete-time 
## stochastic model now used. Deterministic model used to 
## assess R0 and NGM validity.

cd test2
./testscript.sh
cd ..

## same as test2 but using custom mover function
## (note you have to swap the order of the events in Rcpp
## due to the mover function being applied in a different
## order to the standard progression terms in MW)

cd test3
./testscript.sh
cd ..

## same as test3 but with FULL PATHWAYS through stages:
## consequently amended discrete-time stochastic model
## and same R0 leads to different nu

cd test4
./testscript.sh
cd ..

## same as test4 but with a change of contact matrix
## partway through the outbreak---starts off as no
## mixing between ages and then flips to coMix
## deterministic model removed here and no final sizes
## since difficult to calculate when contact matrices 
## change

cd test5
./testscript.sh
cd ..

## same as test4 but with different probabilities of
## transition through different pathways - mainly used
## to test R0 and check mover code

cd test6
./testscript.sh
cd ..

## same as test1 but turning off transmission and using
## different seeds at different time points to
## assess seeding function - this only prints output
## as text file to check manually

cd test7
./testscript.sh
cd ..

## same as test7 but running on full ward structure with
## probabilistic seeding in some wards

cd test8
./testscript.sh
cd ..

## same as test8 but with probabilistic age-seeding

cd test9
./testscript.sh
cd ..

