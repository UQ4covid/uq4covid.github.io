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

## run seeding code
R CMD BATCH --no-restore --slave --no-save R_tools/simulateSeedStates.R

metawards -d model_code/ncov.json -m 2011to2019Data\
    -D model_code/demographics.json\
    --mixer model_code/mix_pathways\
    --mover model_code/move_pathways\
    --input inputs/disease.dat\
    --iterator model_code/iterator\
    -u inputs/user_inputs.txt -o raw_outputs --force-overwrite-output \
    --extractor model_code/ward_extractor\
    --start-date 2020/03/06 --theme simple --nsteps 15 --no-spinner \
    --nthreads 4
