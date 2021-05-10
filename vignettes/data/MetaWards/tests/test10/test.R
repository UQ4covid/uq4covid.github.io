library(tidyverse)
library(patchwork)

## source reconstruct function
source("../../R_tools/dataTools.R")

## loop over age classes
rec <- list()
for(k in 1:3) {
    rec[[k]] <- list()
    for(i in 1:8) {
        ## set up list for runs
        rec[[k]][[i]] <- list()
        ## loop over repeats
        for(j in 1:10) {
            ## establish connection
            system(paste0("bzip2 -dkf raw_outputs", k, "/testx", str_pad(j, 3, pad = "0"), "/age", i, ".db.bz2"))
            con <- DBI::dbConnect(RSQLite::SQLite(), paste0("raw_outputs", k, "/testx", str_pad(j, 3, pad = "0"), "/age", i, ".db"))

            ## extract data
            compact <- tbl(con, "compact") %>%
                collect()
            
            DBI::dbDisconnect(con)
            
            ## reconstruct counts from incidence
            rec[[k]][[i]][[j]] <- select(compact, -ward) %>%
                group_by(day) %>%
                summarise_all(sum) %>%
                arrange(day) %>%
                {reconstruct(
                    compact$Einc, compact$Pinc, compact$I1inc, compact$I2inc, compact$RIinc, 
                    compact$DIinc, compact$Ainc, compact$RAinc,
                    compact$Hinc, compact$RHinc, compact$DHinc
                )} %>%
                magrittr::set_colnames(c(
                    "Einc", "E", "Pinc", "P", "I1inc", "I1", "I2inc", "I2", 
                    "RI", "DI",
                    "Ainc", "A", "RA",
                    "Hinc", "H", "RH", "DH"
                )) %>%
                as_tibble()
            rec[[k]][[i]][[j]]$day <- compact$day
        }
        
        ## collapse to mean and CIs
        rec[[k]][[i]] <- bind_rows(rec[[k]][[i]], .id = "rep") %>%
            select(-ends_with("inc")) %>%
            complete(rep, day = 1:max(.$day)) %>%
            group_by(rep) %>%
            mutate_at(-c(1, 2), ~ifelse(day == 1 & is.na(.), 0, .)) %>%
            fill(names(.)) %>%
            ungroup() %>%
            select(-rep) %>%
            gather(stage, count, -day) %>%
            group_by(day, stage) %>%
            nest() %>%
            mutate(data = map(data, ~{
                data.frame(
                    mn = mean(.$count), 
                    LCI = quantile(.$count, probs = 0.025), 
                    UCI = quantile(.$count, probs = 0.975)
                )
            })) %>%
            unnest(cols = data) %>%
            gather(estimate, value, -day, -stage) %>%
            spread(stage, value) %>%
            select(day, estimate, E, P, I1, DI)
    }
}

## collapse over experiments
rec <- map(rec, ~bind_rows(., .id = "age")) %>%
    bind_rows(.id = "experiment") %>%
    mutate(experiment = paste("experiment", experiment))

## plots
pdf("test.pdf", width = 10, height = 5)
for(i in 1:8) { 
    p <- filter(rec, age == i) %>%
        select(-age) %>%
        gather(stage, value, -day, -experiment, -estimate) %>%
        spread(estimate, value) %>%
        ggplot(aes(x = day)) +
            geom_ribbon(aes(ymin = LCI, ymax = UCI, fill = stage), alpha = 0.3) +
            geom_line(aes(y = mn, colour = stage)) +
            facet_wrap(~experiment) +
            ggtitle(paste0("age", i))
    print(p)
}
dev.off()

