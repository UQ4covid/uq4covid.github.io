#!/bin/bash

mkdir -p raw_outputs

R CMD BATCH --no-restore --slave --no-save convertDesign.R 

export METAWARDSDATA=$HOME/Documents/covid/MetaWardsData

metawards --nproc 24 --nthreads 1 -d model_code/ncov.json -D model_code/demographics.json --mixer model_code/mix_pathways --mover model_code/move_pathways --input inputs/disease.dat -a ExtraSeedsLondon.dat -u model_code/lockdown_states.txt -o raw_outputs --force-overwrite-output --iterator model_code/iterate --extractor model_code/ward_extractor --start-date 2020/01/01 --theme simple --nsteps 100


