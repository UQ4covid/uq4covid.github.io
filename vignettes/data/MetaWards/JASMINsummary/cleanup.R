## R script to query database
library(tidyverse)
library(magrittr)

## set directory to save outputs to and lookup path
filedir <- readLines("filedir.txt")
id <- readLines("id.txt")

## read in job lookup
hashes <- readLines("job_lookup.txt")

## read in output files
output <- map(hashes, function(hash, filedir, id) {
    system(paste0("rm ", filedir, "raw_outputs/", hash, "/output_", id, ".rds"))
}, filedir = filedir, id = id)

## remove output files and uncompressed db
system("rm *.err *.out *.Rout")
system(paste0("rm ", filedir, "raw_outputs/summaries_", id, ".db"))

print("Finished")

