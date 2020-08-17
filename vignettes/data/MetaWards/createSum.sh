#!/bin/bash

## unzip archive
bzip2 -dkf raw_outputs/$1/stages.db.bz2

## run SQL script to produce summary table
sqlite3 "raw_outputs/$1/stages.db" ".read extractSQL.sql"

sqlite3 "raw_outputs/$1/stages.db" ".mode csv" ".headers on" ".output raw_outputs/$1/weeksums.csv" "select * from weeksums;" ".exit"

#mv $1stages.db /gws/nopw/j04/covid19/public/$1stages.db
mkdir -p public/raw_outputs/$1
mv raw_outputs/$1/stages.db public/raw_outputs/$1/
mv raw_outputs/$1/weeksums.csv public/raw_outputs/$1/

