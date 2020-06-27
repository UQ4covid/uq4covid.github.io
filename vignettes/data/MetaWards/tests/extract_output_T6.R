## R script to query database
# library(dbplyr)
# library(RSQLite)
library(tidyverse)

## establish connection
system("bzip2 -dkf raw_outputs/stages.db.bz2")
con <- DBI::dbConnect(RSQLite::SQLite(), "raw_outputs/stages.db")

## extract data
compact <- tbl(con, "compact")

## the commands above don't pull down from the database
## until you run a 'collect()' pull, but it allows
## you to use 'tidyverse'-esque notation, which is converts
## to SQL if you're not used to SQL

## see e.g. https://db.rstudio.com/dplyr

## check all demographics other than genpop are empty
stopifnot(
    select(compact, -day, -ward, -Einc, -E, -Iinc, -I, -RI, -DI) %>%
    collect() %>%
    sum() == 0
)

## plot infection counts
p <- select(compact, day, ward, I) %>%
    group_by(day) %>%
    summarise(I = sum(I)) %>%
    collect() %>%
    ggplot(aes(x = day, y = I)) +
        geom_line() +
        xlab("Day") + ylab("Infections") +
        ggtitle("R0 = 2 Inc./Inf. period = 2 days")
ggsave("infections_T6.pdf", p)

## calculate number of wards infected
wards <- select(compact, ward) %>%
    collect() %>%
    pluck("ward") %>%
    unique()
print(paste0("Number of infected wards = ", length(wards)))

## extract cumulative infections at day 100
genpop <- select(compact, day, ward, Iinc) %>%
    filter(day <= 100) %>%
    arrange(ward, day) %>%
    group_by(ward) %>%
    summarise(Icum = sum(Iinc)) %>%
    collect()

## optimise proportion    
fun <- function(p) abs(1 - exp(-2 * p) - p)
opt <- optimise(fun, interval = c(0, 1))

## extract national level
print(paste0("runs: ", sum(genpop$Icum) / 56082077, "pred: ", opt$minimum))

## remove db
system("rm raw_outputs/stages.db")
