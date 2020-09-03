## R script to query database
library(readr)
library(purrr)
library(dplyr)
library(parallel)

## set path to files and week/output to extract
filedir <- "https://gws-access.jasmin.ac.uk/public/covid19/"
week <- 12
output <- "Hprev"

## read in inputs
design <- read_csv(paste0(filedir, "inputs/design.csv"))
parRanges <- read_csv(paste0(filedir, "inputs/parRanges.csv"))

## run through hashes and extract runs
out <- mclapply(design$output, function(hash, filedir, week, output) {
    out <- read_csv(paste0(path, "raw_outputs/", hash, "/", output, ".csv"))
    if(!is.na(week)) {
        out <- filter(out, week = week) %>% select(-week)
    }
    if(!is.na(output)) {
        out <- select(out, !!output)
    }
    out
}, filedir = filedir, week = week, output = output, mc.cores = ncores)

## bind rows
names(out) <- design$output
out <- bind_rows(out, .id = output)


