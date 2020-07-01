## R script to query database
library(dplyr)
library(purrr)

## may also need to install RSQLite library
## e.g. install.packages("RSQLite")

## read in inputs
design <- readRDS("inputs/design.rds")
parRanges <- readRDS("inputs/parRanges.rds")

## set up path to SQL database
path <- "raw_outputs/uberStages.db"

## establish connection
con <- DBI::dbConnect(RSQLite::SQLite(), path)

## Hprev in the compact table contains the mean hospital counts
## currently week is stored as a character
compact <- tbl(con, "compact")
Hprev_12 <- filter(compact, week == "12") %>%
    select(output, replicate, ward, Hprev) %>%
    collect()

## disconnect from DB
DBI::dbDisconnect(con)

## condense to mean across replicates
Hprev_12 <- group_by(Hprev_12, output, replicate) %>%
    summarise(Hprev = sum(Hprev)) %>%
    summarise(Hprev = mean(Hprev))

