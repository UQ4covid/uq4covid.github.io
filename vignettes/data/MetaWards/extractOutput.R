## load jaspy to get SQLite3 3.26.0 and 
## later version of R (3.5.1)
system("module load jaspy")

## R script to query database
library(dplyr)
library(purrr)
library(stringr)
## later versions of tidyr won't compile
#library(versions)
#install.versions("tidyr", version = "1.0.0")
library(tidyr)
library(parallel)

## may also need to install RSQLite library
## e.g. install.packages("RSQLite")

## read in inputs
design <- readRDS("inputs/design.rds")
parRanges <- readRDS("inputs/parRanges.rds")

## extract data
system.time(output <- map2(design$output, design$repeats, function(hash, reps) {
  
    ## generate output hashes
    hashes <- map2(hash, 1:reps, c)
    
    ## print progress
    print(paste0("Currently evaluating: ", hash))
    
    ## run through hashes and extract runs
    output <- mclapply(hashes, function(hash) {
    
        ## extract replicate and hash
        replicate <- hash[2]
        hash <- hash[1]
        
        ## create path
        path <- paste0(hash, ifelse(replicate > 1, paste0("x", str_pad(replicate, 3, pad = "0")), ""))
        path <- paste0("raw_outputs/", path, "/stages.db")

        ## establish connection
        con <- DBI::dbConnect(RSQLite::SQLite(), path)

        ## Hprev in the weeksums table contains the mean hospital counts
        weeksums <- tbl(con, "weeksums")
        out <- collect(weeksums)

        ## disconnect from DB
        DBI::dbDisconnect(con)
        
        ## return output
        out        
    }, mc.cores = 20)
    
    ## bind rows
    output <- bind_rows(output) %>%
        group_by(ward) %>%
        nest() %>%
        ungroup()
        
    output$data <- mclapply(output$data, function(output) {
        ## extract quantiles (this is faster than using 
        ## multiple "quantile()" calls in "summarise()")
        group_by(output, week) %>%
        summarise_all(~{
            q = list(quantile(., probs = c(0.025, 0.25, 0.5, 0.75, 0.975)))
        }, .groups = "drop") %>%
        mutate_at(-1, function(q) {
            map(q, ~as_tibble(t(.))) %>%
            map(~set_names(., c("q025", "q25", "q5", "q75", "q975")))
        }) %>%
        gather(out, value, -week) %>%
        unnest(cols = value)
    }, mc.cores = 20)
    
    print(paste0(hash, ": Finished"))
    
    ## return outputs
    unnest(output, cols = data) %>%
        mutate(output = hash)
}))

## bind rows
output <- bind_rows(output)

## extract outputs of interest
Hprev <- filter(output, out == "Hprev") %>%
    select(-out)
saveRDS(Hprev, "Hprev.rds")

Cprev <- filter(output, out == "Cprev") %>%
    select(-out)
saveRDS(Cprev, "Cprev.rds")

Deaths <- filter(output, out == "Deaths") %>%
    select(-out)
saveRDS(Deaths, "Deaths.rds")


