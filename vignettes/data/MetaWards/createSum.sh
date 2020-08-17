#!/bin/bash

## unzip archive
bzip2 -dkf $1.bz2

## run SQL script to produce summary table
sqlite3 "$1" ".read extractSQL.sql"

