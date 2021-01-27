#!/bin/bash

mkdir -p raw_outputs

export METAWARDSDATA=$HOME/Documents/covid/MetaWardsData

nprocessors=1
nthreads=4


## test for no progression from R to D
cp diseaseTest.dat disease.dat

## .nu .pE .pEA 
printf '0.2939086 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0 0 0 0 0 0 0 0 ' >> disease.dat
## .pA 
printf '0 0 0 0 0 0 0 0 ' >> disease.dat
## .pI1 .pI1H .pI1I2 
printf '0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ' >> disease.dat
## .pI2
printf '0 0 0 0 0 0 0 0 ' >> disease.dat
## .pH .pHR 
printf '0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ' >> disease.dat
## .lock_1_restrict .lock_2_release .GP_A .GP_H repeats 
printf '0 0 0 0 1' >> disease.dat

## run MetaWards
metawards --nproc $nprocessors --nthreads $nthreads -m 2011to2019Data\
    -d ../model_code/ncov_overall.json\
    -D ../model_code/demographics.json --mixer ../model_code/mix_pathways\
    --mover ../model_code/move_pathways --input disease.dat\
    -u user_inputs.txt -o raw_outputs --force-overwrite-output \
    --iterator ../model_code/iterator --extractor ../model_code/ward_extractor\
    --start-date 2020/01/01 --theme simple --nsteps 10

## run R check
R CMD BATCH --no-save --no-restore --slave extract_output_T1.R



## turn hospitals on
cp diseaseTest.dat disease.dat

## .nu .pE .pEA 
printf '0.2939086 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0 0 0 0 0 0 0 0 ' >> disease.dat
## .pA 
printf '0 0 0 0 0 0 0 0 ' >> disease.dat
## .pI1 .pI1H .pI1I2 
printf '0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0 0 0 0 0 0 0 0 ' >> disease.dat
## .pI2
printf '0 0 0 0 0 0 0 0 ' >> disease.dat
## .pH .pHR 
printf '0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0 0 0 0 0 0 0 0 ' >> disease.dat
## .lock_1_restrict .lock_2_release .GP_A .GP_H repeats 
printf '0 0 0 0 1' >> disease.dat

## run MetaWards
metawards --nproc $nprocessors --nthreads $nthreads -m 2011to2019Data\
    -d ../model_code/ncov_overall.json\
    -D ../model_code/demographics.json --mixer ../model_code/mix_pathways\
    --mover ../model_code/move_pathways --input disease.dat\
    -u user_inputs.txt -o raw_outputs --force-overwrite-output\
    --iterator ../model_code/iterator --extractor ../model_code/ward_extractor\
    --start-date 2020/01/01 --theme simple --nsteps 10

## run R check
R CMD BATCH --no-save --no-restore --slave extract_output_T2.R



## turn asymptomatics on
cp diseaseTest.dat disease.dat

## .nu .pE .pEA 
printf '0.2939086 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 1 1 1 1 1 1 1 1 ' >> disease.dat
## .pA 
printf '0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 ' >> disease.dat
## .pI1 .pI1H .pI1I2 
printf '0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ' >> disease.dat
## .pI2
printf '0 0 0 0 0 0 0 0 ' >> disease.dat
## .pH .pHR 
printf '0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ' >> disease.dat
## .lock_1_restrict .lock_2_release .GP_A .GP_H repeats 
printf '0 0 0 0 1' >> disease.dat

## run MetaWards
metawards --nproc $nprocessors --nthreads $nthreads -m 2011to2019Data\
    -d ../model_code/ncov_overall.json\
    -D ../model_code/demographics.json --mixer ../model_code/mix_pathways\
    --mover ../model_code/move_pathways --input disease.dat\
    -u user_inputs.txt -o raw_outputs --force-overwrite-output\
    --iterator ../model_code/iterator --extractor ../model_code/ward_extractor\
    --start-date 2020/01/01 --theme simple --nsteps 10

## run R check
R CMD BATCH --no-save --no-restore --slave extract_output_T3.R



## turn everything on and turn up R0 to get infections in
## different age-classes
cp diseaseTest.dat disease.dat

## .nu .pE .pEA 
printf '0.4 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 ' >> disease.dat
## .pA 
printf '0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 ' >> disease.dat
## .pI1 .pI1H .pI1I2 
printf '0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.7 0.7 0.7 0.7 0.7 0.7 0.7 0.7 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 ' >> disease.dat
## .pI2
printf '0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 ' >> disease.dat
## .pH .pHR 
printf '0.8 0.8 0.8 0.8 0.8 0.8 0.8 0.8 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 ' >> disease.dat
## .lock_1_restrict .lock_2_release .GP_A .GP_H repeats 
printf '0 0 0.2 0.2 1' >> disease.dat

## run MetaWards
metawards --nproc $nprocessors --nthreads $nthreads -m 2011to2019Data\
    -d ../model_code/ncov_overall.json\
    -D ../model_code/demographics.json --mixer ../model_code/mix_pathways\
    --mover ../model_code/move_pathways --input disease.dat\
    -u user_inputs.txt -o raw_outputs --force-overwrite-output\
    --iterator ../model_code/iterator --extractor ../model_code/ward_extractor\
    --start-date 2020/01/01 --theme simple --nsteps 50

## run R check
R CMD BATCH --no-save --no-restore --slave extract_output_T4.R



## check for saturation across wards with high-ish R0
## just SEIR for simplicity
cp diseaseTest.dat disease.dat

## .nu .pE .pEA 
printf '0.5143401 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0 0 0 0 0 0 0 0 ' >> disease.dat
## .pA 
printf '0 0 0 0 0 0 0 0 ' >> disease.dat
## .pI1 .pI1H .pI1I2 
printf '0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 ' >> disease.dat
## .pI2
printf '0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 ' >> disease.dat
## .pH .pHR 
printf '0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ' >> disease.dat
## .lock_1_restrict .lock_2_release .GP_A .GP_H repeats 
printf '0 0 0 0 1' >> disease.dat

## run MetaWards
metawards --nproc $nprocessors --nthreads $nthreads -m 2011to2019Data\
    -d ../model_code/ncov_overall.json\
    -D ../model_code/demographics.json --mixer ../model_code/mix_pathways\
    --mover ../model_code/move_pathways --input disease.dat\
    -u user_inputs.txt -o raw_outputs --force-overwrite-output\
    --iterator ../model_code/iterator --extractor ../model_code/ward_extractor\
    --start-date 2020/01/01 --theme simple --nsteps 100

## run R check
R CMD BATCH --no-save --no-restore --slave extract_output_T5.R



## check for saturation across wards with high-ish R0
cp diseaseTest.dat disease.dat

## .nu .pE .pEA 
printf '0.5143401 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 ' >> disease.dat
## .pA 
printf '0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 ' >> disease.dat
## .pI1 .pI1H .pI1I2 
printf '0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.25 0.25 0.25 0.25 0.25 0.25 0.25 0.25 ' >> disease.dat
## .pI2
printf '0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 ' >> disease.dat
## .pH .pHR 
printf '0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 ' >> disease.dat
## .lock_1_restrict .lock_2_release .GP_A .GP_H repeats 
printf '0 0 0 0 1' >> disease.dat

## run MetaWards
metawards --nproc $nprocessors --nthreads $nthreads -m 2011to2019Data\
    -d ../model_code/ncov_overall.json\
    -D ../model_code/demographics.json --mixer ../model_code/mix_pathways\
    --mover ../model_code/move_pathways --input disease.dat\
    -u user_inputs.txt -o raw_outputs --force-overwrite-output\
    --iterator ../model_code/iterator --extractor ../model_code/ward_extractor\
    --start-date 2020/01/01 --theme simple --nsteps 100

## run R check
R CMD BATCH --no-save --no-restore --slave extract_output_T6.R


