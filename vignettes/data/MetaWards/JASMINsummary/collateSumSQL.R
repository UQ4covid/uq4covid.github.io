## R script to query database
library(tidyverse)
library(magrittr)

## set directory to save outputs to and lookup path
filedir <- readLines("filedir.txt")
id <- readLines("id.txt")

## read in job lookup
hashes <- readLines("job_lookup.txt")

## set up connection to database
con <- DBI::dbConnect(RSQLite::SQLite(), paste0(filedir, "raw_outputs/summaries_", id, ".db"))

## read in output files
output <- map(hashes, function(hash, filedir, id, con) {
    output <- readRDS(paste0(filedir, "raw_outputs/", hash, "/output_", id, ".rds")) %>%
        mutate(output = hash)
    DBI::dbWriteTable(con, "compact", output, append = TRUE)
}, filedir = filedir, id = id, con = con)

## disconnect from database
DBI::dbDisconnect(con)

## compress file
system(paste0("bzip2 -zf ", filedir, "raw_outputs/summaries_", id, ".db"))

print("Finished")

