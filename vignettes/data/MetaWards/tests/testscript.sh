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
