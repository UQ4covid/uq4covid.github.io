## R script to query database
library(dplyr)
library(purrr)
library(stringr)
library(lubridate)
library(tidyr)
library(parallel)

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

## extract data
output <- map2(design$output[1], design$repeats[1], function(hash, reps) {
  
    ## generate output hashes
    hashes <- map(1:reps, function(i, hash) {
        c(paste0(hash, ifelse(i > 1, paste0("x", str_pad(i, 3, pad = "0")), "")), i)
    }, hash = hash)
    
    ## print progress
    print(paste0("Currently evaluating: ", hash))
    
    ## run through hashes and extract runs
    output <- lapply(hashes, function(hash) {
    
    browser()
    
        ## extract replicate and hash
        replicate <- hash[2]
        path <- hash[1]
        hash <- str_replace(path, "x[:digit:]+$", "")

        ## create path
        path <- paste0("raw_outputs/", path, "/stages.db")
        
        ## run SQL
        cmd <- paste0("sqlite3 \"", path, "\" \".param set :id ", hash, "\" \".param set :rep ", replicate, "\" \".read extractSQL.sql\"")
        system(cmd)
        
    })#, mc.cores = 10)
       
    print(paste0(hashes, ": Finished"))
    
    ## return hash
    hash
})

print("All done.")

