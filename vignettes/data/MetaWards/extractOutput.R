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
        out <- filter(weeksums, week == 12) %>%
            select(ward, Hprev) %>%
            collect()

        ## disconnect from DB
        DBI::dbDisconnect(con)
        
        ## return output
        out        
    }, mc.cores = 20)
    
    ## bind rows
    output <- bind_rows(output, .id = "replicate")
        
    ## extract quantiles (this is faster than using 
    ## multiple "quantile()" calls in "summarise()")
    output <- group_by(output, ward) %>%
        summarise(
            q = list(quantile(Hprev, probs = c(0.025, 0.25, 0.5, 0.75, 0.975)))
        ) %>%
        mutate(q = map(q, ~as_tibble(t(.)))) %>%
        mutate(q = map(q, ~set_names(., c("q025", "q25", "q5", "q75", "q975")))) %>%
        unnest(cols = q) %>%
        mutate(output = hash)
    
    print(paste0(hash, ": Finished"))
    
    ## return outputs
    output
}))

## bind rows
output <- bind_rows(output)

## save output
saveRDS(output, "Hprev_12.rds")


