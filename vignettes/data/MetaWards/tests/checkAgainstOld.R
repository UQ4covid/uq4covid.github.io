## R script to query database
# library(dbplyr)
# library(RSQLite)
library(tidyverse)

## IMPORTANT: this compares outputs from old extractor, commit: 3193748edc9
## with new extractor, commit: 0e64768af
##
## run "convertDesign.R" from both commits, but save the old extractor
## in folder "raw_outputs_origextractor" and the new extractor in
## "raw_outputs". Second time use:
## 
## metawards -c raw_outputs_origextractor/config.yaml
## 
## to rerun with exact same seeds etc. The code below then aims to 
## reconstruct the original outputs from the new outputs to check 
## for validity

## set max time
maxT <- 30

## extract runs
files <- list.files("../raw_outputs/")
files <- files[grep("Ens000", files)]

## extract runs
filesO <- list.files("../raw_outputs_origextractor/")
filesO <- filesO[grep("Ens000", filesO)]

stopifnot(identical(files, filesO))

if(length(files) == 0) {
    files <- "."
}

## source reconstruct function
source("reconstruct.R")

## loop over outputs
for(file in files) {

    ## loop over age classes, unzip and extract day ranges for each ward
    minMaxT <- list(NULL)
    for(i in 1:8) {
        ## unzip databases
        system(paste0("bzip2 -dkf ../raw_outputs/", file, "/stages", i, ".db.bz2"))
        system(paste0("bzip2 -dkf ../raw_outputs_origextractor/", file, "/stages", i, ".db.bz2"))
        
        ## extract data
        con <- DBI::dbConnect(RSQLite::SQLite(), paste0("../raw_outputs/", file, "/stages", i, ".db"))
        compact <- tbl(con, "compact")
        
        ## extract max time point to compare
        minMaxT[[i]] <- select(compact, ward, day) %>%
            distinct() %>%
            collect() %>%
            group_by(ward) %>%
            summarise(min = min(day), max = maxT)
        
        ## disconnect from database
        DBI::dbDisconnect(con)
    }
    ## old extractor only updated stages after a ward was initialised
    ## so mirror that here
    minMaxT <- bind_rows(minMaxT, .id = "stage") %>%
        mutate(stage = as.numeric(stage)) %>%
        complete(stage = 1:8, nesting(ward)) %>%
        group_by(ward) %>%
        mutate(minstage = min(which(min == min(min, na.rm = TRUE)))) %>%
        mutate(min = ifelse(min >= min(min, na.rm = TRUE) & stage < minstage, NA, min)) %>%
        select(-minstage) %>%
        fill(min, .direction = "down") %>%
        mutate(min = ifelse(!is.na(min), min(min, na.rm = TRUE), min(min, na.rm = TRUE) + 1)) %>%
        mutate(max = maxT)

    ## loop over age classes
    for(i in 1:8) {
        ## establish connection
        con <- DBI::dbConnect(RSQLite::SQLite(), paste0("../raw_outputs/", file, "/stages", i, ".db"))
        compact <- tbl(con, "compact")
        
        ## establish connection
        conO <- DBI::dbConnect(RSQLite::SQLite(), paste0("../raw_outputs_origextractor/", file, "/stages", i, ".db"))
        compactO <- tbl(conO, "compact")

        ## the commands above don't pull down from the database
        ## until you run a 'collect()' pull, but it allows
        ## you to use 'tidyverse'-esque notation, which it converts
        ## to SQL if you're not used to SQL

        ## see e.g. https://db.rstudio.com/dplyr
        
        ## extract min/max bounds
        temp <- filter(minMaxT, stage == i)
        
        ## extract and reconstruct
        compact <- collect(compact) %>%
            complete(day = 1:maxT, ward = temp$ward) %>%
            left_join(temp, by = "ward") %>%
            filter(day >= min & day <= max) %>%
            select(-min, -max) %>%
            mutate_all(~ifelse(is.na(.), 0, .)) %>%
            arrange(ward, day)
            
        compact <- group_by(compact, ward) %>%
            nest() %>%
            mutate(data = map(data, ~{
                rec <- reconstruct(
                    .$Einc, .$Iinc, .$Rinc, .$Dinc, .$IAinc, .$RAinc,
                    .$IHinc, .$RHinc, .$DHinc,.$ICinc, .$RCinc, .$DCinc
                ) %>%
                magrittr::set_colnames(c(
                    "Einc", "E", "Iinc", "I", "R", "D",
                    "IAinc", "IA", "RA",
                    "IHinc", "IH", "RH", "DH",
                    "ICinc", "IC", "RC", "DC"
                )) %>%
                as_tibble()
                rec$day <- .$day
                rec
            })) %>%
            unnest(cols = data) %>%
            select(day, everything())
        
        ## extract original
        compactO <- collect(compactO) %>%
            arrange(ward, day)
        
        ## compare reconstruction to original
        stopifnot(all((compact - compactO) == 0))
        
        ## disconnect from databases
        DBI::dbDisconnect(con)
        DBI::dbDisconnect(conO)

        print(paste0("All OK - file ", file, " stage ", i))
    }
    
    ## cleanup
    for(i in 1:8) {
        system(paste0("rm ../raw_outputs/", file, "/stages", i, ".db"))
        system(paste0("rm ../raw_outputs_origextractor/", file, "/stages", i, ".db"))
    }
}
