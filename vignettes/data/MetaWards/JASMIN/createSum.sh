#!/bin/bash

## read in jobs
readarray -t jobs < "job_lookup.txt"

## extract path name
jobname=${jobs[$1-1]}
echo $jobname

## set directory to copy files to
filedir="/gws/nopw/j04/covid19/public/wave0/"

## run R script
cmd="R CMD BATCH --no-restore --no-save --slave '--args ${filedir} ${jobname}' postProc.R postProc_${jobname}.Rout"
eval $cmd

