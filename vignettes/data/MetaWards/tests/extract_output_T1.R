## R script to query database
# library(dbplyr)
# library(RSQLite)
library(tidyverse)

## establish connection
system("bzip2 -dkf raw_outputs/stages3.db.bz2")
con <- DBI::dbConnect(RSQLite::SQLite(), "raw_outputs/stages3.db")

## extract data
compact <- tbl(con, "compact")

## the commands above don't pull down from the database
## until you run a 'collect()' pull, but it allows
## you to use 'tidyverse'-esque notation, which is converts
## to SQL if you're not used to SQL

## see e.g. https://db.rstudio.com/dplyr

## now extract and reconstruct by filling in gaps
compact <- collect(compact)

genpop <- select(compact, day, ward, starts_with("genpop_")) %>%
    complete(day = 1:10, nesting(ward)) %>%
    arrange(ward, day) %>%
    mutate_at(vars(starts_with("genpop_")), ~ifelse(day == 1 & is.na(.), 0, .)) %>%
    fill(starts_with("genpop_"), .direction = "down") %>%
    set_names(c("day", "ward", paste0("stage_", 0:5)))
asymp <- select(compact, day, ward, starts_with("asymp_")) %>%
    gather(stage, value, -day, -ward) %>%
    complete(stage = paste0("asymp_", 0:5), nesting(ward, day), fill =list(value = 0)) %>%
    spread(stage, value) %>%
    complete(day = 1:10, nesting(ward)) %>%
    arrange(ward, day) %>%
    mutate_at(vars(starts_with("asymp_")), ~ifelse(day == 1 & is.na(.), 0, .)) %>%
    fill(starts_with("asymp_"), .direction = "down") %>%
    set_names(c("day", "ward", paste0("stage_", 0:5)))
hospital <- select(compact, day, ward, starts_with("hospital_")) %>%
    gather(stage, value, -day, -ward) %>%
    complete(stage = paste0("hospital_", 0:5), nesting(ward, day), fill =list(value = 0)) %>%
    spread(stage, value) %>%
    complete(day = 1:10, nesting(ward)) %>%
    arrange(ward, day) %>%
    mutate_at(vars(starts_with("hospital_")), ~ifelse(day == 1 & is.na(.), 0, .)) %>%
    fill(starts_with("hospital_"), .direction = "down") %>%
    set_names(c("day", "ward", paste0("stage_", 0:5)))
critical <- select(compact, day, ward, starts_with("critical_")) %>%
    gather(stage, value, -day, -ward) %>%
    complete(stage = paste0("critical_", 0:5), nesting(ward, day), fill =list(value = 0)) %>%
    spread(stage, value) %>%
    complete(day = 1:10, nesting(ward)) %>%
    arrange(ward, day) %>%
    mutate_at(vars(starts_with("critical_")), ~ifelse(day == 1 & is.na(.), 0, .)) %>%
    fill(starts_with("critical_"), .direction = "down") %>%
    set_names(c("day", "ward", paste0("stage_", 0:5)))
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
