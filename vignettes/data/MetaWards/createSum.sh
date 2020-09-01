#!/bin/bash

## read in jobs
readarray -t jobs < "job_lookup.txt"

## extract path name
jobname=${jobs[$1-1]}
echo $jobname

## set directory to copy files to
filedir="/gws/nopw/j04/covid19/public/raw_outputs"

## check directory exists
if [ ! -d "${filedir}/${jobname}" ]; then

    ## unzip archive
    bzip2 -dkf raw_outputs/${jobname}/stages.db.bz2

    sleep 60s

    ## run SQL script to produce summary table
    sqlite3 "raw_outputs/${jobname}/stages.db" ".read extractSQL.sql"

    sleep 60s

    sqlite3 "raw_outputs/${jobname}/stages.db" -cmd ".mode csv" -cmd ".headers on" -cmd ".output raw_outputs/${jobname}/weeksums.csv" "select * from weeksums;" -cmd ".exit"

    sleep 60s

    mkdir -p ${filedir}/${jobname}
    mv raw_outputs/${jobname}/stages.db ${filedir}/${jobname}/
    mv raw_outputs/${jobname}/weeksums.csv ${filedir}/${jobname}/
 else
    echo "${jobfile} directory already exists"
 fi
echo "Complete"

