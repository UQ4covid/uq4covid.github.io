#!/bin/bash

## read in jobs
readarray -t jobs < "job_lookup.txt"

## set up new file
if [ -f "job_check.txt" ]; then
    rm job_check.txt
fi
touch job_check.txt

## extract number of entries of job file
N=${#jobs[@]}

k=1
while [ $k -le $N ]; do
    ## extract path name
    jobname=${jobs[$k-1]}
    
    ## check output files
    jobfile="$1_$k.out"
    completed=$( tail -n 1 ${jobfile} )
    if [ ${completed} == "Complete" ]; then
        jobfile_err="$1_$k.err"
        err=$( cat ${jobfile_err} | wc -l )
        if [ ${err} != "1" ]; then
            echo $k,${jobname} >> job_check.txt
        fi
    else
        echo $k,${jobname} >> job_check.txt
    fi
    ((k++))
done

