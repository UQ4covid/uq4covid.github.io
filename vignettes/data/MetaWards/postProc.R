## R script to query database
library(dplyr)
library(tidyr)
library(readr)

## may also need to install RSQLite library
## e.g. install.packages("RSQLite")

## extract command line arguments
args <- commandArgs(TRUE)
if(length(args) > 0) {
    stopifnot(length(args) == 2)
    filedir <- args[1]
    hash <- args[2]
} else {
    stop("No arguments")
}

## read in ward and week lookups
weeks <- read_csv("week_lookup.csv", col_names = FALSE)
wards <- read_csv("ward_lookup.csv", col_names = FALSE)
colnames(weeks) <- c("day", "week")
colnames(wards) <- c("ward", "week")

## here we are going to extract cumulative hospital counts
## hospital deaths, critical care counts and critical care
## deaths at each week since lockdown

## print progress
print(paste0("Currently evaluating: ", hash))

## create path
path <- paste0("raw_outputs/", hash, "/stages.db")

## unzip DB
system(paste0("bzip2 -dkf ", path, ".bz2"))

## open database connection
con <- DBI::dbConnect(RSQLite::SQLite(), path)

## collect database
output <- tbl(con, "compact") %>%
    collect()
    
## join with week lookup and generate summary measures
output <- inner_join(output, weeks, by = "day") %>%
    group_by(ward, week) %>%
    summarise(Hprev = sum(H) / 7, Cprev = sum(C) / 7, Deaths = max(DH) + max(DC)) %>%
    ungroup()
    
## expand to empty wards
output <- left_join(wards, output, by = c("ward", "week")) %>%
    mutate_all(~replace_na(., 0))

## write table out
system(paste0("mkdir -p ", filedir, "raw_outputs/", hash))
system(paste0("mv ", path, " ", filedir, "raw_outputs/", hash))
write_csv(output, paste0(filedir, "raw_outputs/", hash, "/weeksums.csv"))
 
print(paste0(hash, ": Finished"))

