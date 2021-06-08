library(tidyverse)
library(patchwork)

## source reconstruct function
source("../../R_tools/dataTools.R")

## loop over models
rec <- list()
models <- c("raw_outputs", "raw_outputs1")
for(m in 1:length(models)) {
    ## loop over age classes
    rec[[m]] <- list()
    for(j in 1:10) {
        rec[[m]][[j]] <- list()
        for(i in 1:8) {
            ## establish connection
            system(paste0("bzip2 -dkf ", models[m], "/testx", str_pad(j, 3, pad = "0"), "/age", i, ".db.bz2"))
            con <- DBI::dbConnect(RSQLite::SQLite(), paste0(models[m], "/testx", str_pad(j, 3, pad = "0"), "/age", i, ".db"))

            ## extract data
            compact <- tbl(con, "compact") %>%
                collect()
            
            DBI::dbDisconnect(con)
            
            system(paste0("rm ", models[m], "/testx", str_pad(j, 3, pad = "0"), "/age", i, ".db"))
            
            ## reconstruct counts from incidence
            rec[[m]][[j]][[i]] <- select(compact, day, Einc) %>%
                group_by(day) %>%
                summarise(Einc = sum(Einc)) %>%
                ungroup() %>%
                arrange(day) %>%
                mutate(Ecum = cumsum(Einc))
        }
        
        ## collapse to mean and CIs
        rec[[m]][[j]] <- bind_rows(rec[[m]][[j]], .id = "age") %>%
            complete(age, day = 1:max(.$day)) %>%
            group_by(age) %>%
            mutate_at(-c(1, 2), ~ifelse(day == 1 & is.na(.), 0, .)) %>%
            fill(names(.)) %>%
            ungroup()
    }
    rec[[m]] <- bind_rows(rec[[m]], .id = "rep") %>%
        group_by(age, day) %>%
        summarise(LCI = quantile(Ecum, probs = 0.025), UCI = quantile(Ecum, probs = 0.975), Ecum = mean(Ecum)) %>%
        ungroup()
}
names(rec) <- c("pweekend = 0", "pweekend = 1")
rec <- bind_rows(rec, .id = "model")

pdf("test.pdf", width = 10, height = 5)
p <- ggplot(rec, aes(x = day, y = Ecum)) +
    geom_ribbon(aes(ymin = LCI, ymax = UCI, fill = model), alpha = 0.5) +
    geom_line(aes(colour = model)) +
    facet_wrap(~age)
print(p)
dev.off()

