## R script to query database
library(dplyr)
library(purrr)
library(stringr)
library(ggplot2)

## may also need to install RSQLite library
## e.g. install.packages("RSQLite")

## read in inputs
design <- readRDS("inputs/design.rds")
parRanges <- readRDS("inputs/parRanges.rds")

## read in simulations for each input
## this involves unzipping the DB files
## which is done on the command line at the moment
## need a better solution going forwards

## the commands below do not pull down from the database
## until you run the 'collect()' function, but it allows
## you to use 'tidyverse'-esque notation, which it converts
## to SQL behind the scenes if you're not used to SQL
## see e.g. https://cran.r-project.org/web/packages/dbplyr/vignettes/dbplyr.html

## here we are going to extract cumulative hospital counts
## at day 100

## extract data
cH <- map2(design$output, design$repeats, function(hash, reps) {
    
    ## generate output hashes
    hashes <- map_chr(1:reps, function(i, hash) {
        paste0(hash, ifelse(i > 1, paste0("x", str_pad(i, 3, pad = "0")), ""))
    }, hash = hash)
    
    ## run through hashes and extract runs
    output <- map(hashes, function(hash) {
        ## create path
        path <- paste0("raw_outputs/", hash, "/stages.db")
        
        ## unzip DB
        system(paste0("bzip2 -dk ", path, ".bz2"))
        
        ## establish connection
        con <- DBI::dbConnect(RSQLite::SQLite(), path)
        
        ## stage_2 in the hospital demographic contains the new incidence
        hospital <- tbl(con, "hospital_totals")
        hosp_db <- filter(hospital, day <= 100) %>%
            select(ward, stage_2) %>%
            group_by(ward) %>%
            summarise(cH = sum(stage_2))
        ## collect outcome of query
        hosp <- collect(hosp_db)
        
        ## disconnect from DB
        DBI::dbDisconnect(con)
        
        ## remove DB
        ## using gio trash just in case anything goes wrong
        system(paste0("gio trash ", path))
        
        ## return counts by ward
        hosp
    })
    names(output) <- hashes
    output <- bind_rows(output, .id = "output")
    
    ## remove replicate indexes
    output <- mutate(output, replicate = as.numeric(factor(output))) %>%
        mutate(output = str_replace(output, "x[:digit:]+$", ""))
    stopifnot(length(unique(output$output)) == 1)
    print(paste0(hashes, ": Finished"))
    output
})
cH <- bind_rows(cH)

## now you can join the design file to the outputs
## via the 'output' indicator

## for brevity let's sum hospital cases over all wards
## and then join to the design file for plotting
sims <- group_by(cH, output, replicate) %>%
    summarise(cH = sum(cH)) %>%
    ungroup() %>%
    inner_join(design, by = "output") 

## plot R0 against sims
ggplot(sims, aes(x = r_zero, y = cH, colour = output)) +
    geom_point()
