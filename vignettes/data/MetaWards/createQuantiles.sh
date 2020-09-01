#!/bin/bash

## read in jobs
readarray -t jobs < "job_lookup_R.txt"

## extract path name
jobname=${jobs[$1-1]}
echo $jobname

## set directory to read files from and store results
filedir="/gws/nopw/j04/covid19/public/raw_outputs"

## run R script
cmd="R CMD BATCH --no-restore --no-save --slave '--args ${filedir} ${jobname} $2 $3' extractOutput.R extractOutput_${jobname}.Rout"
eval $cmd

