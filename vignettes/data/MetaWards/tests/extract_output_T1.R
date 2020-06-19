## R script to query database
# library(dbplyr)
# library(RSQLite)
library(tidyverse)

## establish connection
system("bzip2 -d raw_outputs/stages.db.bz2")
con <- DBI::dbConnect(RSQLite::SQLite(), "raw_outputs/stages.db")

## extract data
genpop <- tbl(con, "genpop_totals")
asymp <- tbl(con, "asymp_totals")
hospital <- tbl(con, "hospital_totals")
critical <- tbl(con, "critical_totals")

## the commands above don't pull down from the database
## until you run a 'collect()' pull, but it allows
## you to use 'tidyverse'-esque notation, which is converts
## to SQL if you're not used to SQL

## see e.g. https://db.rstudio.com/dplyr

genpop <- collect(genpop) %>%
    arrange(ward, day)
asymp <- collect(asymp) %>%
    arrange(ward, day)
hospital <- collect(hospital) %>%
    arrange(ward, day)
critical <- collect(critical) %>%
    arrange(ward, day)
results <- list(GP = genpop, A = asymp, H = hospital, C = critical) %>%
    bind_rows(.id = "demo")

## check that there are no removals
stopifnot(all(results$stage_4 == 0))

## check that there are no individuals not in genpop
stopifnot(
    filter(results, demo != "GP") %>%
    select(starts_with("stage")) %>% 
    summarise_all(~all(. == 0)) %>%
    unlist() %>%
    all()
)

## check for entries in death category
stopifnot(any(genpop$stage_5 > 0))

## check for non-decreasing D
stopifnot(
    group_by(genpop, ward) %>%
    nest() %>%
    mutate(data = map_lgl(data, function(x) {
        mutate(x, lg = lag(stage_5)) %>%
        mutate(diff = stage_5 - lg) %>%
        pluck("diff") %>% 
        {all(.[!is.na(.)] >=0)}
    })) %>%
    pluck("data") %>%
    all()
)

filter(genpop, ward == 255)

print("All OK.")
