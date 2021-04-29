## load jaspy to get SQLite3 3.30.1 and 
## later version of R (3.6.3) you will have
## to do this outside of the script
# system("module load jaspy")
print("HAVE YOU LOADED jaspy?")

## set name of directory to append to save outputs to and startdate
## (directory appended to "/gws/nopw/j04/covid19/public")
filedir <- "wave0"
startdate <- "09/02/2020"
ndays <- 41

## check filedir and reformat if necessary
filedir <- gsub("^/", "", filedir)
filedir <- paste0("/gws/nopw/j04/covid19/public/", filedir)
filedir <- gsub("//$", "/", paste0(filedir, "/"))

## make public directory if it's not there
system(paste0("mkdir -p ", filedir))
system(paste0("chmod -R u=rwx,g=rwx,o=r ", filedir))

## write to external file
system(paste0("echo ", filedir, " > filedir.txt"))

## load libraries
library(tidyverse)
library(lubridate)

## create weeks from startdate
startdate <- dmy(startdate)
dates <- startdate + 0:ndays
tweeks <- week(dates)
WEEKS <- unique(tweeks)

## create week/day lookup table
week_lookup <- data.frame(day = as.numeric(dates - startdate), date = dates, week = tweeks)
week_lookup <- filter(week_lookup, week %in% WEEKS)
write.table(week_lookup, "week_lookup.csv", row.names = FALSE, col.names = FALSE, sep = ",")

## create ward lookup table
ward_lookup <- expand.grid(ward = 1:8071, week = unique(week_lookup$week))
write.table(ward_lookup, "ward_lookup.csv", row.names = FALSE, col.names = FALSE, sep = ",")

## read in inputs
inputs <- readRDS("../inputs/inputs.rds")
parRanges <- readRDS("../inputs/parRanges.rds")

## write out as csvs
system(paste0("mkdir -p ", filedir, "inputs"))
write_csv(inputs, paste0(filedir, "inputs/inputs.csv"))
write_csv(parRanges, paste0(filedir, "inputs/parRanges.csv"))

## this next part creates the script file to pass to SLURM

## extract data
paths <- map2(inputs$output, inputs$repeats, function(hash, reps) {
  
    ## generate output hashes
    hashes <- map2(hash, 1:reps, c)
    
    ## print progress
    print(paste0("Currently evaluating: ", hash))
    
    ## run through hashes and extract runs
    output <- map_chr(hashes, function(hash, append) {
    
        ## extract replicate and hash
        replicate <- hash[2]
        hash <- hash[1]

        ## create path
        path <- paste0(hash, ifelse(append, paste0("x", str_pad(replicate, 3, pad = "0")), ""))
        path
    }, append = (reps > 1))
    
    print(paste0(hash, ": Finished"))
    
    ## return hash
    output
})
paths <- reduce(paths, c)

## write number of jobs to file
code <- readLines("submit_job_template.sbatch")
code <- gsub("RANGES", paste0("1-", length(paths)), code)
code <- gsub("FILEDIR", filedir, code)
writeLines(code, "submit_job.sbatch")

## write csv to query
paths <- data.frame(path = paths)
write.table(paths, "job_lookup.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)

print("All done.")

