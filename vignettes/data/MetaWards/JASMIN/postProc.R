## R script to query database
library(tidyverse)
library(magrittr)

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

## source in dataTools (for reconstruct function)
source("../R_tools/dataTools.R")

## read in ward and week lookups
weeks <- read_csv("week_lookup.csv", col_names = FALSE)
wards <- read_csv("ward_lookup.csv", col_names = FALSE)
colnames(weeks) <- c("day", "date", "week")
colnames(wards) <- c("ward", "week")
weeks <- mutate(weeks, date = as.Date(date, format = "%Y-%m-%d"))

## here we are going to extract cumulative hospital counts
## hospital deaths, critical care counts and critical care
## deaths at each week since lockdown

## print progress
print(paste0("Currently evaluating: ", hash))

## extract files
files <- list.files(paste0("../raw_outputs/", hash))
files <- files[grep(glob2rx("age*.db.bz2"), files)]

for(file in files) {

    ## set path
    path <- paste0("../raw_outputs/", hash, "/", file)
    path <- gsub(".bz2", "", path)

    ## unzip DB
    system(paste0("bzip2 -dkf ", path, ".bz2"))

    ## open database connection
    con <- DBI::dbConnect(RSQLite::SQLite(), path)

    ## collect database
    compact <- tbl(con, "compact") %>%
        collect()
        
    ## disconnect from database
    DBI::dbDisconnect(con)
        
    ## reconstruct counts from incidence
    output <- arrange(compact, day) %>%
        group_by(ward) %>%
        nest() %>%
        mutate(data = map(data, ~{
            output <- reconstruct(
                .$Einc, .$Pinc, .$I1inc, .$I2inc, .$RIinc, 
                .$DIinc, .$Ainc, .$RAinc,
                .$Hinc, .$RHinc, .$DHinc
            ) %>%
            magrittr::set_colnames(c(
                "Einc", "E", "Pinc", "P", "I1inc", "I1", "I2inc", "I2", 
                "RI", "DI",
                "Ainc", "A", "RA",
                "Hinc", "H", "RH", "DH"
            )) %>%
            as_tibble()
            output$day <- .$day
            output
        })) %>%
        unnest(cols = data) %>%
        select(day, ward, everything())

    ## join with week lookup and generate summary measures
    output <- inner_join(output, weeks, by = "day") %>%
        group_by(ward, week) %>%
        summarise(Hprev = sum(H) / 7, Hdeaths = max(DH), Cdeaths = max(DI)) %>%
        ungroup()
        
    ## expand to empty wards
    output <- left_join(wards, output, by = c("ward", "week")) %>%
        mutate_all(~replace_na(., 0)) %>%
        left_join(
            group_by(weeks, week) %>%
            summarise(date = max(date))
        , by = "week")
        
    ## age class
    agec <- strsplit(file, "\\.")[[1]][1]

    ## write table out
    system(paste0("mkdir -p ", filedir, "raw_outputs/", hash))
    system(paste0("mv ", path, " ", filedir, "raw_outputs/", hash))
    write_csv(output, paste0(filedir, "raw_outputs/", hash, "/weeksums_", agec, ".csv"))
    
    print(paste0("Done: ", agec))
}
 
print(paste0(hash, ": Finished"))

