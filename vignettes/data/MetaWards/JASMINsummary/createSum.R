## R script to query database
library(tidyverse)
library(magrittr)

## may also need to install RSQLite library
## e.g. install.packages("RSQLite")

## extract command line arguments
args <- commandArgs(TRUE)
if(length(args) > 0) {
    stopifnot(length(args) == 3)
    filedir <- args[1]
    hash <- args[2]
    id <- args[3]
} else {
    stop("No arguments")
}

## source in dataTools (for reconstruct function)
source("../R_tools/dataTools.R")

## Add trust and wards for DW
trusts <- read_csv("../inputs/trust19Lookup.csv")
ward19 <- read_csv("../inputs/Ward19_Lookup.csv")
TrustLookup <- inner_join(trusts, ward19, by = c("code" = "WD19CD")) %>%
    rename(ward = FID) %>%
    select(ward, trustId)

## read in ward and week lookups
weeks <- read_csv("../JASMINsetup/week_lookup.csv", col_names = FALSE)
wards <- read_csv("../JASMINsetup/ward_lookup.csv", col_names = FALSE)
colnames(weeks) <- c("day", "date", "week")
colnames(wards) <- c("ward", "week")
weeks <- mutate(weeks, date = as.Date(date, format = "%Y-%m-%d"))

## combine lookups
lookup <- full_join(weeks, wards, by = "week") %>%
    select(-date)

## here we are going to extract cumulative hospital counts
## hospital deaths, critical care counts and critical care
## deaths at each week since lockdown

## print progress
print(paste0("Currently evaluating: ", hash))

## extract files
files <- list.files(paste0(filedir, "raw_outputs/", hash))
files <- files[grep(glob2rx("age*.db"), files)]

## create output file
output <- map(files, function(file, filedir, hash, lookup) {

    ## set up path
    path <- paste0(filedir, "raw_outputs/", hash, "/", file)

    ## open database connection
    con <- DBI::dbConnect(RSQLite::SQLite(), path)

    ## collect databases
    compact <- tbl(con, "compact") %>%
        collect()
    compact_ini <- tbl(con, "compact_ini") %>%
        collect() %>%
        set_names(colnames(compact))
        
    ## disconnect from database
    DBI::dbDisconnect(con)
    
    ## combine databases together and complete 
    ## initial conditions
    compact <- rbind(compact_ini, compact) %>%
        complete(day = 1, ward)
    compact[is.na(compact)] <- 0
        
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
    ## expanding data sets where necessary        
    output <- select(output, day, ward, Hinc, H, DH, DI) %>%
        full_join(lookup, by = c("day", "ward")) %>%
        arrange(ward, day) %>%
        mutate(Hinc = ifelse(is.na(Hinc), 0, Hinc)) %>%
        mutate(H = ifelse(day == 0 & is.na(H), 0, H)) %>%
        mutate(DH = ifelse(day == 0 & is.na(DH), 0, DH)) %>%
        mutate(DI = ifelse(day == 0 & is.na(DI), 0, DI)) %>%
        group_by(ward) %>%
        fill(H) %>%
        fill(DH) %>%
        fill(DI) %>%
        mutate(cumH = cumsum(Hinc)) %>%
        group_by(ward, week) %>%
        summarise(cumH = max(cumH), Hprev_mn = sum(H) / 7, cumHD = max(DH), cumCD = max(DI)) %>%
        ungroup()

       ## collate to trust
       output <- inner_join(output, TrustLookup, by = "ward") %>%
            group_by(week, trustId) %>%
            summarise(cumH = sum(cumH), Hprev_mn = sum(Hprev_mn), cumHD = sum(cumHD), cumCD = sum(cumCD)) %>%
            ungroup()
       print(head(output)) 
    
    ## return summaries
    output
}, filedir = filedir, hash = hash, lookup = lookup)

## create uber file
names(output) <- files
output <- bind_rows(output, .id = "age") %>%
    mutate(age = gsub(".db", "", age)) %>%
    mutate(age = gsub("age", "", age)) %>%
    mutate(age = as.numeric(age))

## write table out
saveRDS(output, paste0(filedir, "raw_outputs/", hash, "/output_", id, ".rds"))
 
print(paste0(hash, ": Finished"))

