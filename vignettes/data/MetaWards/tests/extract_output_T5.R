## R script to query database
# library(dbplyr)
# library(RSQLite)
library(tidyverse)

## source reconstruct function
source("reconstruct.R")

## set up proportions in each age class
ageprop <- c(0.06, 0.154, 0.154, 0.134, 0.128, 0.134, 0.105, 0.131)

## loop over age classes
for(i in 1:8) {
    ## establish connection
    system(paste0("bzip2 -dkf raw_outputs/stages", i, ".db.bz2"))
    con <- DBI::dbConnect(RSQLite::SQLite(), paste0("raw_outputs/stages", i, ".db"))

    ## extract data
    compact <- tbl(con, "compact") %>%
        collect()
    
    DBI::dbDisconnect(con)
    
    ## reconstruct counts from incidence
    compact <- group_by(compact, ward) %>%
        nest() %>%
        mutate(data = map(data, ~{
            rec <- reconstruct(
                .$Einc, .$Iinc, .$Rinc, .$Dinc, .$IAinc, .$RAinc,
                .$IHinc, .$RHinc, .$DHinc,.$ICinc, .$RCinc, .$DCinc
            ) %>%
                magrittr::set_colnames(c(
                    "Einc", "E", "Iinc", "I", "R", "D",
                    "IAinc", "IA", "RA",
                    "IHinc", "IH", "RH", "DH",
                    "ICinc", "IC", "RC", "DC"
                )) %>%
                as_tibble()
            rec$day <- .$day
            rec
        })) %>%
        unnest(cols = data) %>%
        select(day, everything()) %>%
        ungroup()
    
    ## check all demographics other than genpop are empty
    stopifnot(
        select(compact, -day, -ward, -Einc, -E, -Iinc, -I, -R, -D) %>%
        sum() == 0
    )

    ## plot infection counts
    p <- select(compact, day, ward, I) %>%
        group_by(day) %>%
        summarise(I = sum(I)) %>%
        ggplot(aes(x = day, y = I)) +
            geom_line() +
            xlab("Day") + ylab("Infections") +
            ggtitle("R0 = 3.5 Inc./Inf. period = 2 days")
    ggsave(paste0("infections_T5_stage", i, ".pdf"), p)

    ## calculate number of wards infected
    wards <- select(compact, ward) %>%
        pluck("ward") %>%
        unique()
    print(paste0("Number of infected wards = ", length(wards)))

    ## extract cumulative infections at day 100
    genpop <- select(compact, day, ward, Iinc) %>%
        filter(day <= 100) %>%
        arrange(ward, day) %>%
        group_by(ward) %>%
        summarise(Icum = sum(Iinc))

    ## optimise proportion    
    fun <- function(p) abs(1 - exp(-3.5 * p) - p)
    opt <- optimise(fun, interval = c(0, 1))

    ## extract national level
    print(paste0("runs stage ", i, ": ", sum(genpop$Icum) / (56082077 * ageprop[i]), " pred: ", opt$minimum))

    ## remove db
    system(paste0("rm raw_outputs/stages", i, ".db"))
}
