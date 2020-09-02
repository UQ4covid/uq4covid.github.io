#!/bin/bash

## read in jobs
readarray -t jobs < "job_lookup.txt"

## extract path name
jobname=${jobs[$1-1]}
echo $jobname

## set directory to copy files to
filedir="/gws/nopw/j04/covid19/public/raw_outputs"

## copy files to directory
cp week_lookup.csv raw_outputs/${jobname}
cp ward_lookup.csv raw_outputs/${jobname}
cp extractSQL.sql raw_outputs/${jobname}
cp readCSV.sql raw_outputs/${jobname}

## change directory 
cd raw_outputs/${jobname}

## remove old file
if [ -f "stages.db" ]; then
    rm stages.db
fi
if [ -f "stages_${jobname}.db" ]; then
    rm stages_${jobname}.db
fi

## unzip archive
bzip2 -dk stages.db.bz2
mv stages.db stages_${jobname}.db

## run SQL scripts to produce summary table
sqlite3 "stages_${jobname}.db" -cmd ".read readCSV.sql" &&
sqlite3 "stages_${jobname}.db" -cmd ".read extractSQL.sql" &&
sqlite3 "stages_${jobname}.db" -cmd ".mode csv" -cmd ".headers on" -cmd ".output weeksums.csv" "select * from weeksums;" -cmd ".exit"

## move to public repo
mkdir -p ${filedir}/${jobname}
mv stages_${jobname}.db ${filedir}/${jobname}/
mv weeksums.csv ${filedir}/${jobname}/

## remove extraneous files
rm week_lookup.csv
rm ward_lookup.csv
rm extractSQL.sql
rm readCSV.sql

mv ${filedir}/${jobname}/stages_${jobname}.db ${filedir}/${jobname}/stages.db

#echo "Complete"

