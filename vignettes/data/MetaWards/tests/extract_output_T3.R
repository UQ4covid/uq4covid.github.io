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

## check that there are no individuals in any
## class they shouldn't be
stopifnot(
    filter(results, demo != "A") %>%
    select(stage_2, stage_3, stage_4, stage_5) %>% 
    summarise_all(~all(. == 0)) %>%
    unlist() %>%
    all()
)
stopifnot(
    filter(results, demo != "GP") %>%
    select(stage_0, stage_1) %>% 
    summarise_all(~all(. == 0)) %>%
    unlist() %>%
    all()
)

## check that there are infections
stopifnot(any(asymp$stage_2 > 0))
stopifnot(any(asymp$stage_3 > 0))

## check that there are no deaths
stopifnot(all(asymp$stage_5 == 0))

## check for entries in removal category
stopifnot(any(asymp$stage_4 > 0))

## check for non-decreasing R
stopifnot(
    group_by(asymp, ward) %>%
    nest() %>%
    mutate(data = map_lgl(data, function(x) {
        mutate(x, lg = lag(stage_4)) %>%
        mutate(diff = stage_4 - lg) %>%
        pluck("diff") %>% 
        {all(.[!is.na(.)] >=0)}
    })) %>%
    pluck("data") %>%
    all()
)

print("genpop")
filter(genpop, ward == 255)

print("asymp")
filter(asymp, ward == 255)

print("All OK.")
