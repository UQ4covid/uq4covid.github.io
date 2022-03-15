## load jasr reminder
# system("module load jasr")
print("HAVE YOU LOADED jasr?")

## set name of directory to save outputs
wave <- 1

## read in input file
pars <- readRDS(paste0("wave", wave, "/disease.rds"))

## write number of jobs to file
code <- readLines("submit_job_template.sbatch")
code <- gsub("RANGES", paste0("1-", nrow(pars), code)
code <- gsub("FILEDIR", paste0("wave", wave), code)
writeLines(code, "submit_job.sbatch")

## write csv to query
write.table(data.frame(job = 1:nrow(pars)), "job_lookup.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)

print("All done.")

