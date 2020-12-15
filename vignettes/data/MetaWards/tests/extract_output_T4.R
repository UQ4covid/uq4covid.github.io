## R script to query database
# library(dbplyr)
# library(RSQLite)
library(tidyverse)

## source reconstruct function
source("reconstruct.R")

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
            select(day, everything()) %>%
            ungroup()
        
        ## now extract and reconstruct by filling in gaps
        compact <- rename(compact, genpop_0 = Einc, genpop_1 = E, genpop_2 = Iinc, genpop_3 = I, genpop_4 = R, genpop_5 = D) %>%
            rename(asymp_2 = IAinc, asymp_3 = IA, asymp_4 = RA) %>%
            rename(hospital_2 = IHinc, hospital_3 = IH, hospital_4 = RH, hospital_5 = DH) %>%
            rename(critical_2 = ICinc, critical_3 = IC, critical_4 = RC, critical_5 = DC)
    
        genpop <- select(compact, day, ward, starts_with("genpop_")) %>%
            complete(day = 1:10, nesting(ward)) %>%
            arrange(ward, day) %>%
            mutate_at(vars(starts_with("genpop_")), ~ifelse(day == 1 & is.na(.), 0, .)) %>%
            fill(starts_with("genpop_"), .direction = "down") %>%
            set_names(c("day", "ward", paste0("stage_", 0:5)))
        asymp <- select(compact, day, ward, starts_with("asymp_")) %>%
            gather(stage, value, -day, -ward) %>%
            complete(stage = paste0("asymp_", 0:5), nesting(ward, day), fill =list(value = 0)) %>%
            spread(stage, value) %>%
            complete(day = 1:10, nesting(ward)) %>%
            arrange(ward, day) %>%
            mutate_at(vars(starts_with("asymp_")), ~ifelse(day == 1 & is.na(.), 0, .)) %>%
            fill(starts_with("asymp_"), .direction = "down") %>%
            set_names(c("day", "ward", paste0("stage_", 0:5)))
        hospital <- select(compact, day, ward, starts_with("hospital_")) %>%
            gather(stage, value, -day, -ward) %>%
            complete(stage = paste0("hospital_", 0:5), nesting(ward, day), fill =list(value = 0)) %>%
            spread(stage, value) %>%
            complete(day = 1:10, nesting(ward)) %>%
            arrange(ward, day) %>%
            mutate_at(vars(starts_with("hospital_")), ~ifelse(day == 1 & is.na(.), 0, .)) %>%
            fill(starts_with("hospital_"), .direction = "down") %>%
            set_names(c("day", "ward", paste0("stage_", 0:5)))
        critical <- select(compact, day, ward, starts_with("critical_")) %>%
            gather(stage, value, -day, -ward) %>%
            complete(stage = paste0("critical_", 0:5), nesting(ward, day), fill =list(value = 0)) %>%
            spread(stage, value) %>%
            complete(day = 1:10, nesting(ward)) %>%
            arrange(ward, day) %>%
            mutate_at(vars(starts_with("critical_")), ~ifelse(day == 1 & is.na(.), 0, .)) %>%
            fill(starts_with("critical_"), .direction = "down") %>%
            set_names(c("day", "ward", paste0("stage_", 0:5)))
        results <- list(GP = genpop, A = asymp, H = hospital, C = critical) %>%
            bind_rows(.id = "demo")
    
        ## check that there are no individuals in any
        ## class they shouldn't be
        stopifnot(
            filter(results, demo != "GP") %>%
            select(stage_0, stage_1) %>% 
            summarise_all(~all(. == 0)) %>%
            unlist() %>%
            all()
        )
    
        ## check that there are infections, removals and deaths
        ## except in asymptomatics
        stopifnot(
            filter(results, demo != "A") %>%
            select(demo, stage_2, stage_3, stage_4, stage_5) %>%
            group_by(demo) %>%
            summarise_all(~any(. > 0)) %>%
            gather(stage, value, -demo) %>%
            pluck("value") %>%
            all()
        )
        ## check that there are infections and removals
        ## in asymptomatics
        stopifnot(
            filter(results, demo == "A") %>%
            select(stage_2, stage_3, stage_4) %>%
            summarise_all(~any(. > 0)) %>%
            gather(stage, value) %>%
            pluck("value") %>%
            all()
        )
        ## check for no deaths in asymptomatics
        stopifnot(
            filter(results, demo == "A") %>%
            select(stage_5) %>%
            summarise_all(~all(. == 0)) %>%
            pluck("stage_5")
        )
        
        ## check for non-decreasing R
        stopifnot(
            mutate(results, lg = lag(stage_4)) %>%
            mutate(diff = stage_4 - lg) %>%
            mutate(diff = ifelse(day == 1, NA, diff)) %>%
            pluck("diff") %>% 
                {all(.[!is.na(.)] >= 0)}
        )
        
        ## check for non-decreasing D
        stopifnot(
            mutate(results, lg = lag(stage_5)) %>%
            mutate(diff = stage_5 - lg) %>%
            mutate(diff = ifelse(day == 1, NA, diff)) %>%
            pluck("diff") %>% 
            {all(.[!is.na(.)] >= 0)}
        )
    
        print("genpop")
        filter(genpop, ward == 255)
    
        print("asymp")
        filter(asymp, ward == 255)
    
        print("hospital")
        filter(hospital, ward == 255)
    
        print("critical")
        filter(critical, ward == 255)
    
        print(paste0("All OK - stage ", i))
    } else {
        print(paste0("No infections - stage ", i))
    }
}
