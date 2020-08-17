## load jaspy to get SQLite3 3.26.0 and 
## later version of R (3.5.1) you will have
## to do this outside of the script
# system("module load jaspy")
print("HAVE YOU LOADED jaspy?")

## load libraries
library(dplyr)
library(purrr)
library(stringr)
## need to install before-CCTZ version of lubridate
## or else it won't compile
## also then need devtools, so had to install
## earlier version of a dependency
#library(versions)
#install.versions("later", version = "0.7.5")
#library(devtools)
#devtools::install_github("tidyverse/lubridate@before-CCTZ")
library(lubridate)

## may also need to install RSQLite library
## e.g. install.packages("RSQLite")

## create weeks from 1st January 2020
startdate <- dmy("01/01/2020")
dates <- startdate + 0:177

lockdownDate1 <- dmy("21/03/2020")
lockdownDate2 <- dmy("13/05/2020")
wld1 <- week(lockdownDate1)
wld2 <- week(lockdownDate2)
tweeks <- week(dates)

## lockdown 
lockdown1 <- dmy("21/03/2020")
lockdown2 <- dmy("13/05/2020")

## extract weeks after lockdown
WEEKS <- unique(tweeks[dates >= lockdown1])
WEEKS <- c(min(WEEKS) - 1, WEEKS)

## create week/day lookup table
week_lookup <- data.frame(day = as.numeric(dates - startdate), week = tweeks)
week_lookup <- filter(week_lookup, week %in% WEEKS)
write.csv(week_lookup, "week_lookup.csv", row.names = FALSE)

## create ward lookup table
ward_lookup <- expand.grid(ward = 1:8588, week = unique(week_lookup$week))
write.csv(ward_lookup, "ward_lookup.csv", row.names = FALSE)

## read in inputs
design <- readRDS("inputs/design.rds")
parRanges <- readRDS("inputs/parRanges.rds")

## here we are going to extract cumulative hospital counts
## hospital deaths, critical care counts and critical care
## deaths at each week since lockdown

## this next part creates the script file to pass to LOTUS

## extract data
paths <- map2(design$output, design$repeats, function(hash, reps) {
  
    ## generate output hashes
    hashes <- map2(hash, 1:reps, c)
    
    ## print progress
    print(paste0("Currently evaluating: ", hash))
    
    ## run through hashes and extract runs
    output <- map_chr(hashes, function(hash) {
    
        ## extract replicate and hash
        replicate <- hash[2]
        hash <- hash[1]

        ## create path
        path <- paste0(hash, ifelse(replicate > 1, paste0("x", str_pad(replicate, 3, pad = "0")), ""))
        path <- paste0("./createSum.sh ", path)
        path
    })
    
    print(paste0(hash, ": Finished"))
    
    ## return hash
    output
})
paths <- reduce(paths, c)

## write paths to file
code <- readLines("submit_job_template.sh")
code <- gsub("REPLACE_STUFF", paste(paths, collapse = "\n"), code)
writeLines(code, "submit_job.bsub")

print("All done.")

