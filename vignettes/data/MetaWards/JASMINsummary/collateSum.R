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
    readRDS(paste0(filedir, "raw_outputs/", hash, "/output_", id, ".rds"))
}, filedir = filedir, id = id)

## create uber file
names(output) <- hashes
output <- bind_rows(output, .id = "output")

## write table out
write_csv(output, paste0(filedir, "raw_outputs/summaries_", id, ".csv"), col_names = TRUE, row_names = FALSE)

## compress file
system(paste0("bzip2 -zf ", filedir, "raw_outputs/summaries_", id, ".csv"))
 
print("Finished")

