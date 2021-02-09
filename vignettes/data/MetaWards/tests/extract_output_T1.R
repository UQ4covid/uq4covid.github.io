## R script to query database
# library(dbplyr)
# library(RSQLite)
library(tidyverse)

## source reconstruct function
source("../R_tools/dataTools.R")

## loop over age classes
for(i in 1:8) {
    ## establish connection
    system(paste0("bzip2 -dkf raw_outputs/age", i, ".db.bz2"))
    con <- DBI::dbConnect(RSQLite::SQLite(), paste0("raw_outputs/age", i, ".db"))

    ## extract data
    compact <- tbl(con, "compact") %>%
        collect()
    
    DBI::dbDisconnect(con)
    
    if(nrow(compact) > 0) {
        ## reconstruct counts from incidence
        compact <- group_by(compact, ward) %>%
            nest() %>%
            mutate(data = map(data, ~{
                rec <- reconstruct(
                    .$Einc, .$Pinc, .$I1inc, .$I2inc, .$RIinc, .$DIinc, .$Ainc, .$RAinc,
                    .$Hinc, .$RHinc, .$DHinc
                ) %>%
                magrittr::set_colnames(c(
                    "Einc", "E", "Pinc", "P", "I1inc", "I1", "I2inc", "I2", 
                    "RI", "DI",
                    "Ainc", "A", "RA",
                    "Hinc", "H", "RH", "DH"
                )) %>%
                as_tibble()
                rec$day <- .$day
                rec
            })) %>%
            unnest(cols = data) %>%
            select(day, everything()) %>%
            ungroup()
        
        ## now extract and reconstruct by filling in gaps    
        genpop <- select(compact, day, ward, Einc, E, Pinc, P, I1inc, I1, I2inc, I2, RI, DI) %>%
            complete(day = 1:10, nesting(ward)) %>%
            arrange(ward, day) %>%
            mutate_at(-c(1:2), ~ifelse(day == 1 & is.na(.), 0, .)) %>%
            fill(-c(1:2), .direction = "down")
        asymp <- select(compact, day, ward, Ainc, A, RA) %>%
            complete(day = 1:10, nesting(ward)) %>%
            arrange(ward, day) %>%
            mutate_at(-c(1:2), ~ifelse(day == 1 & is.na(.), 0, .)) %>%
            fill(-c(1:2), .direction = "down")
        hospital <- select(compact, day, ward, Hinc, H, RH, DH) %>%
            complete(day = 1:10, nesting(ward)) %>%
            arrange(ward, day) %>%
            mutate_at(-c(1:2), ~ifelse(day == 1 & is.na(.), 0, .)) %>%
            fill(-c(1:2), .direction = "down")
    
        ## check that there are no individuals in any
        ## class/demographic that they shouldn't be in
        stopifnot(
            select(asymp, -day, -ward) %>%
            summarise_all(~all(. == 0)) %>%
            unlist() %>%
            all()
        )
        stopifnot(
            select(hospital, -day, -ward) %>%
            summarise_all(~all(. == 0)) %>%
            unlist() %>%
            all()
        )
    
        ## check that there are no removals or I2
        stopifnot(all(genpop$RI == 0))
        stopifnot(all(genpop$I2 == 0))
    
        ## check for entries in death category
        stopifnot(any(genpop$DI >= 0))
    
        ## check for non-decreasing D
        stopifnot(
            group_by(genpop, ward) %>%
            nest() %>%
            mutate(data = map_lgl(data, function(x) {
                mutate(x, lg = lag(DI)) %>%
                mutate(diff = DI - lg) %>%
                pluck("diff") %>% 
                {all(.[!is.na(.)] >=0)}
            })) %>%
            pluck("data") %>%
            all()
        )
    
        ## filter(genpop, ward == 255)
    
        print(paste0("All OK - stage ", i))
    } else {
        print(paste0("No infections - stage ", i))
    }
}
