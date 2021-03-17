#!/bin/bash

## set up outputs folder
mkdir -p raw_outputs

## run design code
R CMD BATCH --no-restore --slave --no-save convertDesign.R 

## set path to MetaWardsData repository
export METAWARDSDATA=$HOME/Documents/covid/MetaWardsData
  
## run model
metawards -d model_code/ncov.json -m 2011to2019Data\
    -D model_code/demographics.json\
    --mixer model_code/mix_pathways\
    --mover model_code/move_pathways\
    --input inputs/disease.dat\
    --iterator model_code/iterator\
    -u inputs/user_inputs.txt -o raw_outputs --force-overwrite-output \
    --extractor model_code/ward_extractor\
    --start-date 2020/02/22 --theme simple --nsteps 28
