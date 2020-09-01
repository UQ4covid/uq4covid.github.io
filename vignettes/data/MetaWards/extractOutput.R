## R script to query database
library(readr)
library(purrr)
library(dplyr)
library(stringr)
library(tidyr)
library(parallel)

## extract command line arguments
args <- commandArgs(TRUE)
if(length(args) > 0) {
    stopifnot(length(args) == 4)
    filedir <- args[1]
    hash <- args[2]
    reps <- as.numeric(args[3])
    ncores <- as.numeric(args[4])
} else {
    stop("No arguments")
}

## generate output hashes
hashes <- map2(hash, 1:reps, c)

## print progress
print(paste0("Currently evaluating: ", hash))

## run through hashes and extract runs
output <- mclapply(hashes, function(hash, filedir) {

    ## extract replicate and hash
    replicate <- hash[2]
    hash <- hash[1]
    
    ## create path
    path <- paste0(hash, ifelse(replicate > 1, paste0("x", str_pad(replicate, 3, pad = "0")), ""))
    path <- paste0(filedir, "/", path, "/weeksums.csv")

    ## Hprev in the weeksums table contains the mean hospital counts
    out <- read_csv(path)
    
    ## return output
    out        
}, filedir = filedir, mc.cores = ncores)

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
}, mc.cores = ncores)

## return outputs
output <- unnest(output, cols = data) %>%
    mutate(output = hash)

## extract outputs of interest
Hprev <- filter(output, out == "Hprev") %>%
    select(-out)
saveRDS(Hprev, paste0(filedir, "/", hash, "_Hprev.rds"))

Cprev <- filter(output, out == "Cprev") %>%
    select(-out)
saveRDS(Cprev, paste0(filedir, "/", hash, "_Cprev.rds"))

Deaths <- filter(output, out == "Deaths") %>%
    select(-out)
saveRDS(Deaths, paste0(filedir, "/", hash, "_Deaths.rds"))

print(paste0(hash, ": Finished"))


