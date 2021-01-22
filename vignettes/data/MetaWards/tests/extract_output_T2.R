## R script to query database
# library(dbplyr)
# library(RSQLite)
library(tidyverse)

## source reconstruct function
source("../R_tools/dataTools.R")

## loop over age classes
for(i in 1:8) {
    ## establish connection
    system(paste0("bzip2 -dkf raw_outputs/stages", i, ".db.bz2"))
    con <- DBI::dbConnect(RSQLite::SQLite(), paste0("raw_outputs/stages", i, ".db"))

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
                    .$Einc, .$I1inc, .$I2inc, .$Rinc, .$Dinc, .$IAinc, .$RAinc,
                    .$IHinc, .$RHinc, .$DHinc
                ) %>%
                magrittr::set_colnames(c(
                    "Einc", "E", "I1inc", "I1", "I2inc", "I2", 
                    "R", "D",
                    "IAinc", "IA", "RA",
                    "IHinc", "IH", "RH", "DH"
                )) %>%
                as_tibble()
                rec$day <- .$day
                rec
            })) %>%
            unnest(cols = data) %>%
            select(day, everything()) %>%
            ungroup()
        
        ## now extract and reconstruct by filling in gaps    
        genpop <- select(compact, day, ward, Einc, E, I1inc, I1, I2inc, I2, R, D) %>%
            complete(day = 1:10, nesting(ward)) %>%
            arrange(ward, day) %>%
            mutate_at(-c(1:2), ~ifelse(day == 1 & is.na(.), 0, .)) %>%
            fill(-c(1:2), .direction = "down")
        asymp <- select(compact, day, ward, IAinc, IA, RA) %>%
            complete(day = 1:10, nesting(ward)) %>%
            arrange(ward, day) %>%
            mutate_at(-c(1:2), ~ifelse(day == 1 & is.na(.), 0, .)) %>%
            fill(-c(1:2), .direction = "down")
        hospital <- select(compact, day, ward, IHinc, IH, RH, DH) %>%
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
    
        ## check that there are no removals or I2
        stopifnot(all(genpop$R == 0))
        stopifnot(all(genpop$I2 == 0))
        stopifnot(all(hospital$RH == 0))
    
        ## check for entries in death category
        stopifnot(any(genpop$D >= 0))
        stopifnot(any(hospital$DH >= 0))
    
        ## check for non-decreasing D
        stopifnot(
            group_by(genpop, ward) %>%
            nest() %>%
            mutate(data = map_lgl(data, function(x) {
                mutate(x, lg = lag(D)) %>%
                mutate(diff = D - lg) %>%
                pluck("diff") %>% 
                {all(.[!is.na(.)] >=0)}
            })) %>%
            pluck("data") %>%
            all()
        )
        stopifnot(
            group_by(hospital, ward) %>%
            nest() %>%
            mutate(data = map_lgl(data, function(x) {
                mutate(x, lg = lag(DH)) %>%
                mutate(diff = DH - lg) %>%
                pluck("diff") %>% 
                {all(.[!is.na(.)] >=0)}
            })) %>%
            pluck("data") %>%
            all()
        )
    
        print("genpop")
        filter(genpop, ward == 255)
    
        print("hospital")
        filter(hospital, ward == 255)
    
        print(paste0("All OK - stage ", i))
    } else {
        print(paste0("No infections - stage ", i))
    }
}
