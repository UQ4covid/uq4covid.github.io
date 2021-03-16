## load libraries
library(tidyverse)
library(patchwork)

## read in data
deaths <- read_csv("trustDeathsAdmissionsIcuExtract20200624.csv", guess_max = 30000) %>%
    filter(statistic == "death") %>%
    filter(date <= "2020-04-01")

## examine cases by date
p1 <- deaths %>%
    group_by(date) %>%
    summarise(deaths = sum(value)) %>%
    ggplot(aes(x = date, y = deaths)) +
        geom_line() + xlab("Date") + ylab("Deaths")
p2 <- deaths %>%
    group_by(date) %>%
    summarise(deaths = sum(value)) %>%
    mutate(deaths = cumsum(deaths)) %>%
    ggplot(aes(x = date, y = deaths)) +
        geom_line() + xlab("Date") + ylab("Cumulative Deaths")
p <- p1 + p2
ggsave("deathCounts.pdf", p)

## extract first 100 deaths and amalgamate to proportion
## of deaths per trust
seeds <- deaths %>%
    group_by(date) %>%
    summarise(deaths = sum(value)) %>%
    mutate(deaths = cumsum(deaths)) %>%
    filter(deaths <= 100) %>%
    select(-deaths) %>%
    inner_join(deaths, by = "date") %>%
    filter(value > 0) %>%
    group_by(code) %>%
    summarise(deaths = sum(value)) %>%
    mutate(prop = deaths / sum(deaths)) %>%
    select(-deaths)

## now read in ward to trust lookup
trust19Lookup <- read_csv("../hospitalCatchments/trust19Lookup.csv")

## check trust names match
anti_join(seeds, trust19Lookup, by = c("code" = "trustId"))

print("SOME TRUSTS DON'T MATCH")

## temporarily remove mismatches for testing
seeds <- semi_join(seeds, trust19Lookup, by = c("code" = "trustId"))
seeds <- mutate(seeds, prop = prop / sum(prop))

## set path to MetaWardsData
path <- paste0("../../../../../../MetaWardsData/model_data/2011to2019Data/")

## load workers and players
workers <- read_delim(paste0(path, "WorkSize19.dat"), col_names = FALSE, delim = " ")
players <- read_delim(paste0(path, "PlaySize19.dat"), col_names = FALSE, delim = " ")
popsize <- full_join(workers, players, by = "X1") %>%
    replace_na(list(X2.x = 0, X2.y = 0)) %>%
    mutate(popsize = X2.x + X2.y) %>%
    select(-X2.x, -X2.y)

## join to Ward19 lookup
Ward19Lookup <- read_csv(paste0(path, "Ward19_Lookup.csv")) %>%
    inner_join(popsize, by = c("FID" = "X1")) %>%
    inner_join(trust19Lookup, by = c("WD19CD" = "code")) %>%
    group_by(trustId) %>%
    mutate(prop = popsize / sum(popsize)) %>%
    ungroup()

## subset by trusts with initial infections
seeds <- inner_join(Ward19Lookup, seeds, by = c("trustId" = "code")) %>%
    mutate(prop = prop.x * prop.y) %>%
    select(FID, prop)
stopifnot(sum(seeds$prop) == 1)

## write seeds file
write_csv(seeds, "../../inputs/ward_seeds.csv", col_names = FALSE)

