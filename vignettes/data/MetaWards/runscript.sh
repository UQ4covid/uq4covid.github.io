#!/bin/bash

mkdir -p raw_outputs

export METAWARDSDATA=$HOME/Documents/covid/MetaWardsData

metawards --nproc 24 --nthreads 1 -d ncov.json -D demographics.json --mixer mix_pathways --mover move_pathways --input inputs/disease.dat -a ExtraSeedsLondon.dat -u lockdown_states.txt -o raw_outputs --force-overwrite-output --iterator iterate --start-date 2020/01/01 --theme simple --nsteps 10


