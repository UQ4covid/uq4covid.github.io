#!/bin/bash

## firstly run standard model with workers and players
./testscript1.sh

## now run model turning all workers to players via custom iterator 
./testscript2.sh

## now aggregate worker and player data
R CMD BATCH --no-restore --no-save --slave mergeData.R

## copy to MetaWardsData
mv test /home/tj/Documents/covid/MetaWardsData/model_data/

## run model with new data but model_code iterator
./testscript3.sh

## remove test data
rm -r /home/tj/Documents/covid/MetaWardsData/model_data/test

## combine results together
R CMD BATCH --no-restore --no-save --slave test.R

