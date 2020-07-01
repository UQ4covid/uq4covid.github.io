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

## extract Sundays from week before lockdown
WEEKS <- unique(tweeks[dates >= lockdown1])

## read in inputs
design <- readRDS("inputs/design.rds")
parRanges <- readRDS("inputs/parRanges.rds")

## here we are going to extract cumulative hospital counts
## hospital deaths, critical care counts and critical care
## deaths at each week since lockdown

## write to external files
dir.create("outputs", showWarnings = FALSE)

## extract data
output <- map2(design$output, design$repeats, function(hash, reps, WEEKS, dates) {
  
    ## generate output hashes
    hashes <- map_chr(1:reps, function(i, hash) {
    paste0(hash, ifelse(i > 1, paste0("x", str_pad(i, 3, pad = "0")), ""))
    }, hash = hash)
    
    print(paste0("Currently evaluating: ", hash))
    
    ## run through hashes and extract runs
    output <- mclapply(hashes, function(hash, WEEKS, dates) {

        ## create path
        path <- paste0("raw_outputs/", hash, "/stages.db")

        ## unzip DB
        system(paste0("bzip2 -dkf ", path, ".bz2"))

        ## establish connection
        con <- DBI::dbConnect(RSQLite::SQLite(), path)

        ## extract counts from "compact" table
        compact <- tbl(con, "compact")
        output <- dplyr::select(compact, day, ward, H, C, DH, DC) %>%
          collect()

        ## Add variable indicating week
        output  <- mutate(output, week = week(dates[day])) %>%
            filter(week %in% WEEKS) %>%
            group_by(ward, week) %>%
            summarise(Hprev = mean(H), Cprev = mean(C), Deaths = max(DC) + max(DH)) %>%
            ungroup() %>%
            complete(ward = 1:8588, week = WEEKS, fill = list(Hprev = 0, Cprev = 0, Deaths = 0))

        ## disconnect from DB
        DBI::dbDisconnect(con)

        ## remove DB
        system(paste0("rm ", path))

        ## return counts
        output
    }, WEEKS = WEEKS, dates = dates, mc.cores = 4)
    names(output) <- hashes
    output <- bind_rows(output, .id = "output")

    ## remove replicate indexes
    output <- mutate(output, replicate = as.numeric(factor(output))) %>%
        mutate(output = str_replace(output, "x[:digit:]+$", ""))
    stopifnot(length(unique(output$output)) == 1)
    
    ## extract individual outputs at individual weeks
    mclapply(WEEKS, function(wk, output) {
        temp <- filter(output, week == wk)
        vars <- c("Hprev", "Cprev", "Deaths")
        temp <- map(vars, function(out, output) {
            select(output, output, ward, replicate, value = matches(out))
        }, output = temp)
        ## save outputs
        for(i in 1:length(vars)) {
            saveRDS(temp[[i]], paste0("outputs/", hash, "_", vars[i], "_", wk, ".rds"))
        }
    }, output = output, mc.cores = 4)
    
    print(paste0(hashes, ": Finished"))
    ## return hash
    hash
}, WEEKS = WEEKS, dates = dates)

## concatenate data
out_names <- expand.grid(c("Hprev", "Cprev", "Deaths"), WEEKS)
output <- map2(out_names[, 1], out_names[, 2], function(output_name, week, hashes) {
    ## loop over design points
    output <- mclapply(hashes, function(hash, output_name, week) {
    
        path <- paste0("outputs/", hash, "_", output_name, "_", week, ".rds")
        print(paste0("Currently evaluating: ", hash))
        
        if(!file.exists(path)) {
            stop(paste0(path, " doesn't exist"))
        } else {    
            output <- readRDS(path)
        }
        print(paste0(path, ": Finished"))
        
        ## return output
        output
    }, output_name = output_name, week = week, mc.cores = 4)
    ## bind rows
    output <- bind_rows(output)
    saveRDS(output, paste0("outputs/", output_name, "_", week, ".rds"))
}, hashes = design$output)

## clean up
output <- map2(out_names[, 1], out_names[, 2], function(output_name, week, hashes) {
    ## remove old files
    mclapply(hashes, function(hash, output_name, week) {
        path <- paste0("outputs/", hash, "_", output_name, "_", week, ".rds")
        system(paste0("rm ", path))
    }, output_name = output_name, week = week, mc.cores = 4)
}, hashes = design$output)

print("All done.")

