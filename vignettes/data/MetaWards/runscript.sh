#!/bin/bash

mkdir -p raw_outputs

export METAWARDSDATA=$HOME/Documents/covid/MetaWardsData

metawards -d ncov -D demographics.json --mixer mix_pathways --move move_pathways --input inputs/disease.dat -a ExtraSeedsLondon.dat -i inputs/move_mix.dat -o raw_outputs --force-overwrite-output --iterator /home/tj/Documents/covid/uq4covid/model_config/metawards/iterate --start-date 2020/01/01 --theme simple


