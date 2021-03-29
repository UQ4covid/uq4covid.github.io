#!/bin/bash

## read in jobs
readarray -t jobs < "job_lookup.txt"

## extract path name
jobname=${jobs[$1-1]}
echo $jobname

## run R script
cmd="R CMD BATCH --no-restore --no-save --slave '--args $2 ${jobname}' unZIP.R unZIP_${jobname}.Rout"
eval $cmd

