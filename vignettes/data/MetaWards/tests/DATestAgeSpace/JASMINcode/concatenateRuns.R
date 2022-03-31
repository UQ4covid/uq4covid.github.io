## load jasr reminder
# system("module load jasr")
print("HAVE YOU LOADED jasr?")

## load libraries
library(tidyverse)

## set name of directory to search for outcomes
wave <- 1

## read in input file
pars <- readRDS(paste0("../wave", wave, "/disease.rds"))

## extract log-likelihoods
ll <- map_dbl(1:nrow(pars), function(i, wave) {
    print(i)
    if(file.exists(paste0("../wave", wave, "/runs_md_", i, ".rds"))) {
        run <- readRDS(paste0("../wave", wave, "/runs_md_", i, ".rds"))$ll
    } else {
        run <- NA
    }
    run
}, wave = wave)

## check all runs have completed
if(any(is.na(ll))) {
    cat("Missing runs:\n")
    print(which(is.na(ll)))
    stop("Stopped")
}
#write.table(data.frame(ind = which(is.na(ll))), "job_lookup.txt", row.names = FALSE, col.names = FALSE)

## cleanup
map(1:nrow(pars), function(i, wave) {
    system(paste0("rm ../wave", wave, "/wave", wave, "Runs_", i, ".Rout"))
}, wave = wave)

## save
saveRDS(ll, paste0("../wave", wave, "/ll.rds"))
