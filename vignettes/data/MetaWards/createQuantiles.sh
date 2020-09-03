#!/bin/bash

## read in jobs
readarray -t jobs < "quantile_lookup.txt"

## extract path name
jobname=${jobs[$1-1]}
echo $jobname

## set directory to read files from and store results
filedir="/gws/nopw/j04/covid19/public/wave0/"

## run R script
cmd="R CMD BATCH --no-restore --no-save --slave '--args ${filedir} ${jobname} $2 $3' extractQuantiles.R extractQuantiles_${jobname}.Rout"
eval $cmd

