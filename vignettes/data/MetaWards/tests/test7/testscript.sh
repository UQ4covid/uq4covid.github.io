#!/bin/bash

## set up path to data
mkdir -p raw_outputs
export METAWARDSDATA=$HOME/Documents/covid/MetaWardsData

## set number of threads and processors
nprocessors=1
nthreads=1

## set up user inputs
cp ../user_inputs_default.txt user_inputs.txt
printf '\n.ward_seed_filename = \"ward_seeds.csv\"\n' >> user_inputs.txt
printf '\n.age_seed_filename = \"age_seeds.csv\"\n' >> user_inputs.txt
printf '\n.ini_states_filename = \"seeds.csv\"\n' >> user_inputs.txt
printf '\n\n.contact_matrix1_filename = \"contact_matrix.csv\"\n' >> user_inputs.txt
printf '\n\n.contact_matrix2_filename = \"contact_matrix.csv\"\n' >> user_inputs.txt

## set up contact matrix
rm contact_matrix.csv
touch contact_matrix.csv
printf '1,0,0,0,0,0,0,0\n' >> contact_matrix.csv
printf '0,1,0,0,0,0,0,0\n' >> contact_matrix.csv
printf '0,0,1,0,0,0,0,0\n' >> contact_matrix.csv
printf '0,0,0,1,0,0,0,0\n' >> contact_matrix.csv
printf '0,0,0,0,1,0,0,0\n' >> contact_matrix.csv
printf '0,0,0,0,0,1,0,0\n' >> contact_matrix.csv
printf '0,0,0,0,0,0,1,0\n' >> contact_matrix.csv
printf '0,0,0,0,0,0,0,1' >> contact_matrix.csv

## set up seeds
rm ward_seeds.csv
touch ward_seeds.csv
printf '1,1,1' >> ward_seeds.csv
rm age_seeds.csv
touch age_seeds.csv
printf '1,1\n' >> age_seeds.csv
printf '2,0\n' >> age_seeds.csv
printf '3,0\n' >> age_seeds.csv
printf '4,0\n' >> age_seeds.csv
printf '5,0\n' >> age_seeds.csv
printf '6,0\n' >> age_seeds.csv
printf '7,0\n' >> age_seeds.csv
printf '8,0\n' >> age_seeds.csv
rm seeds.csv
touch seeds.csv
printf 'Ens0000x001,1,10,0,0,0,0,0,0,0,0,0,0\n' >> seeds.csv
printf 'Ens0000x002,1,5,0,2,0,0,0,0,0,0,0,0\n' >> seeds.csv
printf 'Ens0000x003,1,3,0,0,0,0,5,0,0,0,0,0\n' >> seeds.csv
   
## SEPID model
cp ../diseaseTest.dat disease.dat

# .pE .pEP
printf '0.3 0.3 0.3 0.3 0.3 0.3 0.3 0.3 1 1 1 1 1 1 1 1 ' >> disease.dat
## .pP 
printf '0.3 0.3 0.3 0.3 0.3 0.3 0.3 0.3 ' >> disease.dat
## .pA 
printf '0 0 0 0 0 0 0 0 ' >> disease.dat
## .pI1 .pI1H .pI1D 
printf '0.3 0.3 0.3 0.3 0.3 0.3 0.3 0.3 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 ' >> disease.dat
## .pI2
printf '0 0 0 0 0 0 0 0 ' >> disease.dat
## .pH .pHD 
printf '0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 ' >> disease.dat
## .lock_1_restrict .lock_2_release 
printf '1 1 ' >> disease.dat
## beta[1] beta[2] beta[3] beta[6]
printf '0 0 0 0 ' >> disease.dat
## pweekend repeats output
printf '0 3 Ens0000' >> disease.dat

## run MetaWards
rm nohup.out
nohup metawards --nproc $nprocessors --nthreads $nthreads -m single -P 100000\
    -d ncov_age.json\
    -D demographics_age.json --mixer ../../model_code/mix_pathways\
    --input disease.dat\
    --iterator ../../model_code/iterator\
    -u user_inputs.txt -o raw_outputs --force-overwrite-output \
    --extractor ../../model_code/ward_extractor\
    --start-date 2019/12/31 --theme simple --nsteps 3
    
## run R script
R CMD BATCH --no-restore --no-save --slave test.R 
    
