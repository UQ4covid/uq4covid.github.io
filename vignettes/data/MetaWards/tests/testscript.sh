#!/bin/bash

mkdir -p raw_outputs

export METAWARDSDATA=$HOME/Documents/covid/MetaWardsData



## test for no progression from R to D
cp diseaseTest.dat disease.dat

## beta[2] beta[3] .pE .pEA 
printf '1 1 0.5 0 ' >> disease.dat
## .pA 
printf '0 ' >> disease.dat
## .pI .pIH .pIR 
printf '0.5 0 0 ' >> disease.dat
## .pH .pHC .pHR 
printf '0 0 0 ' >> disease.dat
## .pC .pCR
printf '0 0 ' >> disease.dat
## .lock_1_restrict .lock_2_release .GP_A .GP_H .GP_C repeats 
printf '0 0 0 0 0 1' >> disease.dat

## run MetaWards
metawards --nproc 24 --nthreads 1 -d ../model_code/ncov.json -D ../model_code/demographics.json --mixer ../model_code/mix_pathways --mover ../model_code/move_pathways --input disease.dat -a ExtraSeedsLondon.dat -u ../model_code/lockdown_states.txt -o raw_outputs --force-overwrite-output --iterator ../model_code/iterate --extractor ../model_code/ward_extractor --start-date 2020/01/01 --theme simple --nsteps 10

## run R check
R CMD BATCH --no-save --no-restore --slave extract_output_T1.R



## turn hospitals on
cp diseaseTest.dat disease.dat

## beta[2] beta[3] .pE .pEA 
printf '1 1 0.5 0 ' >> disease.dat
## .pA 
printf '0 ' >> disease.dat
## .pI .pIH .pIR 
printf '0.5 0.5 0 ' >> disease.dat
## .pH .pHC .pHR 
printf '0.5 0 0 ' >> disease.dat
## .pC .pCR
printf '0 0 ' >> disease.dat
## .lock_1_restrict .lock_2_release .GP_A .GP_H .GP_C repeats 
printf '0 0 0 0 0 1' >> disease.dat

## run MetaWards
metawards --nproc 24 --nthreads 1 -d ../model_code/ncov.json -D ../model_code/demographics.json --mixer ../model_code/mix_pathways --mover ../model_code/move_pathways --input disease.dat -a ExtraSeedsLondon.dat -u ../model_code/lockdown_states.txt -o raw_outputs --force-overwrite-output --iterator ../model_code/iterate --extractor ../model_code/ward_extractor --start-date 2020/01/01 --theme simple --nsteps 10

## run R check
R CMD BATCH --no-save --no-restore --slave extract_output_T2.R



## turn asymptomatics on
cp diseaseTest.dat disease.dat

## beta[2] beta[3] .pE .pEA 
printf '1 1 0.5 1 ' >> disease.dat
## .pA 
printf '0.5 ' >> disease.dat
## .pI .pIH .pIR 
printf '0 0 0 ' >> disease.dat
## .pH .pHC .pHR 
printf '0 0 0 ' >> disease.dat
## .pC .pCR
printf '0 0 ' >> disease.dat
## .lock_1_restrict .lock_2_release .GP_A .GP_H .GP_C repeats 
printf '0 0 0 0 0 1' >> disease.dat

## run MetaWards
metawards --nproc 24 --nthreads 1 -d ../model_code/ncov.json -D ../model_code/demographics.json --mixer ../model_code/mix_pathways --mover ../model_code/move_pathways --input disease.dat -a ExtraSeedsLondon.dat -u ../model_code/lockdown_states.txt -o raw_outputs --force-overwrite-output --iterator ../model_code/iterate --extractor ../model_code/ward_extractor --start-date 2020/01/01 --theme simple --nsteps 10

## run R check
R CMD BATCH --no-save --no-restore --slave extract_output_T3.R



## turn everything on
cp diseaseTest.dat disease.dat

## beta[2] beta[3] .pE .pEA 
printf '1 1 0.5 0.5 ' >> disease.dat
## .pA 
printf '0.5 ' >> disease.dat
## .pI .pIH .pIR 
printf '0.5 0.7 0.1 ' >> disease.dat
## .pH .pHC .pHR 
printf '0.8 0.8 0.1 ' >> disease.dat
## .pC .pCR
printf '0.5 0.1 ' >> disease.dat
## .lock_1_restrict .lock_2_release .GP_A .GP_H .GP_C repeats 
printf '0 0 0.2 0.2 0.2 1' >> disease.dat

## run MetaWards
metawards --nproc 24 --nthreads 1 -d ../model_code/ncov.json -D ../model_code/demographics.json --mixer ../model_code/mix_pathways --mover ../model_code/move_pathways --input disease.dat -a ExtraSeedsLondon.dat -u ../model_code/lockdown_states.txt -o raw_outputs --force-overwrite-output --iterator ../model_code/iterate --extractor ../model_code/ward_extractor --start-date 2020/01/01 --theme simple --nsteps 50

## run R check
R CMD BATCH --no-save --no-restore --slave extract_output_T4.R

