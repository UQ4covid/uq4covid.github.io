## load jaspy to get SQLite3 3.30.1 and 
## later version of R (3.6.3) you will have
## to do this outside of the script
# system("module load jaspy")
print("HAVE YOU LOADED jaspy?")

## extract and check unique ID
id <- basename(getwd())
id <- strsplit(id, "JASMINsummary_")[[1]]
if(length(id) != 2) {
    stop("Folder name not of form 'JASMINsummary_*'")
}
id <- id[2]
if(id == "") {
    stop("Folder name not of form 'JASMINsummary_*'")
}

## check folder exists
filedir <- readLines("../JASMINsetup/filedir.txt")
if(!dir.exists(filedir)) {
    stop(paste0(filedir, " directory does not exist; you need to run \"JASMINsetup\" code"))
}

## write to output file
system(paste0("echo ", filedir, " > filedir.txt"))
system(paste0("echo ", id, " > id.txt"))

## load libraries
library(tidyverse)

## read in inputs
inputs <- readRDS("../inputs/inputs.rds")

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

## check these paths existed in setup runs
pathsSetup <- readLines("../JASMINsetup/job_lookup.txt")
if(!all(paths %in% pathsSetup)) {
    stop("Some 'paths' have not been extracted")
}

## write number of jobs to file
code <- readLines("submit_job_template.sbatch")
code <- gsub("RANGES", paste0("1-", length(paths)), code)
code <- gsub("FILEDIR", filedir, code)
code <- gsub("USER", id, code)
writeLines(code, "submit_job.sbatch")

## write csv to query
paths <- data.frame(path = paths)
write.table(paths, "job_lookup.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)

print("All done.")

