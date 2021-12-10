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
        
        ## bind together
        compact <- rbind(compact_ini, compact) %>%
            arrange(day)
        
        ## reconstruct counts from incidence
        rec[[i]][[j]] <- reconstruct(
                    compact$Einc, compact$Pinc, compact$I1inc, compact$I2inc, compact$RIinc, 
                    compact$DIinc, compact$Ainc, compact$RAinc,
                    compact$Hinc, compact$RHinc, compact$DHinc
                ) %>%
                magrittr::set_colnames(c(
                    "Einc", "E", "Pinc", "P", "I1inc", "I1", "I2inc", "I2", 
                    "RI", "DI",
                    "Ainc", "A", "RA",
                    "Hinc", "H", "RH", "DH"
                )) %>%
                as_tibble()
        rec[[i]][[j]]$day <- compact$day
        print(paste("i = ", i, ", j = ", j))
        print(rec[[i]][[j]])
    }
}

