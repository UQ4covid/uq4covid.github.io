#!/bin/bash

## read in jobs
readarray -t jobs < "job_lookup.txt"

## set directory to copy files to
filedir="/gws/nopw/j04/covid19/public/raw_outputs"

## set up new file
if [ -f "job_check.txt" ]; then
    rm job_check.txt
fi
touch job_check.txt

## extract number of entries of job file
N=${#jobs[@]}

k=0
while [ $k -lt $N ]; do
    ## extract path name
    jobname=${jobs[$k]}
    
    ## check directory exists
    if [ ! -d "${filedir}/${jobname}" ]; then
        echo "${jobname}" >> job_check.txt
    fi
    ((k++))
done

