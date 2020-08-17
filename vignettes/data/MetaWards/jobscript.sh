#!/bin/bash
#PBS -l walltime=24:00:00
#PBS -l select=8:ncpus=64:mem=64GB
# The above sets 8 nodes with 64 cores each

# source the skeleton bashrc to make sure we have all variables that we need
source /home/metawards/OFFICIAL_SENSITIVE/bash_template

source $METADIR/envs/metawards-latest/bin/activate

# change into the directory from which this job was submitted
cd $PBS_O_WORKDIR

mkdir -p raw_outputs

#R CMD BATCH --no-restore --slave --no-save convertDesign.R 

metawards -d model_code/ncov.json -D model_code/demographics.json \
          --mixer model_code/mix_pathways --mover model_code/move_pathways \
          --input inputs/disease.dat -a ExtraSeedsLondon.dat \
          -u model_code/lockdown_states.txt -o raw_outputs --force-overwrite-output \
          --iterator model_code/iterate --extractor model_code/ward_extractor \
          --start-date 2020/01/01 --theme simple --nsteps 177 --no-spinner \
          --nthreads 4
