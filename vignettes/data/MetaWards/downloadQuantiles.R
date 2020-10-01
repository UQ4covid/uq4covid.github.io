## R script to query database
library(readr)
library(purrr)
library(dplyr)
library(parallel)

## set path to files and week/output to extract
filedir <- "https://gws-access.jasmin.ac.uk/public/covid19/wave0/"
week <- 12
output <- "Hprev"
ncores <- 24

## read in inputs
design <- read_csv(paste0(filedir, "inputs/design.csv"))
parRanges <- read_csv(paste0(filedir, "inputs/parRanges.csv"))

## run through hashes and extract runs
out <- mclapply(design$output, function(hash, filedir, weekReq, output) {
    out <- read_csv(paste0(filedir, "raw_outputs/", hash, "/", output, ".csv"))
    if(!is.na(weekReq[1])) {
        out <- filter(out, week %in% weekReq)
    }
    out
}, filedir = filedir, weekReq = week, output = output, mc.cores = ncores)

## bind rows
out <- bind_rows(out)


