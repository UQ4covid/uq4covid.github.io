library(tidyverse)
library(patchwork)

## source reconstruct function
source("../../R_tools/dataTools.R")

## loop over age classes
rec <- list()
for(i in 1:8) {
    ## set up list for runs
    rec[[i]] <- list()
    ## loop over repeats
    for(j in 1:3) {
        ## establish connection
        system(paste0("bzip2 -dkf raw_outputs/Ens0000x", str_pad(j, 3, pad = "0"), "/age", i, ".db.bz2"))
        con <- DBI::dbConnect(RSQLite::SQLite(), paste0("raw_outputs/Ens0000x", str_pad(j, 3, pad = "0"), "/age", i, ".db"))

        ## extract data
        compact <- tbl(con, "compact") %>%
            collect()
        compact_ini <- tbl(con, "compact_ini") %>%
            collect() %>%
            set_names(colnames(compact))
            
        DBI::dbDisconnect(con)
            
        if(nrow(compact) > 0 | nrow(compact_ini) > 0) {
            ## remove wards with no events
            colnames(compact_ini) <- colnames(compact)
            compact <- rbind(compact_ini, compact) %>%
                complete(day = 1, ward)
            compact[is.na(compact)] <- 0
            compact <- arrange(compact, ward, day)
            
            ## reconstruct counts from incidence
            rec[[i]][[j]] <- group_by(compact, ward) %>%
                summarise(data = reconstruct(Einc, Pinc, I1inc, I2inc, RIinc, DIinc, Ainc, RAinc, Hinc, RHinc, DHinc) %>%
                magrittr::set_colnames(c(
                    "Einc", "E", "Pinc", "P", "I1inc", "I1", "I2inc", "I2", 
                    "RI", "DI",
                    "Ainc", "A", "RA",
                    "Hinc", "H", "RH", "DH"
                )) %>%
                as_tibble() %>%
                mutate(day = day) %>%
                list()
            ) %>%
            unnest(cols = data) %>%
            select(day, ward, everything()) %>%
            arrange(day, ward) %>%
            {filter(., rowSums(.[, -c(1, 2)]) != 0)} %>%
            as.data.frame()
        } else {
            rec[[i]][[j]] <- "No infection"
        }
    }
}

for(i in 1:3) {
    for(j in 1:8) {
        print(paste("i = ", i, ", j = ", j))
        print(rec[[j]][[i]])
     }
 }

