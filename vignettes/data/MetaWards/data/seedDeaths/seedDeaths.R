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

## now read in ward to trust lookup
trust19Lookup <- read_csv("../hospitalCatchments/trust19Lookup.csv")

## extract deaths before 14th March 
seeds <- deaths %>%
    arrange(date) %>%
    group_by(date, code) %>%
    summarise(deaths = sum(value)) %>%
    filter(date <= "2020-03-14") %>%
    filter(deaths > 0) %>%
    ungroup()
    
## check trust names match
temp <- anti_join(seeds, trust19Lookup, by = c("code" = "trustId"))

print(paste0(sum(temp$deaths), " DEATHS DON'T MATCH TO TRUSTS IN LOOKUP"))

## temporarily remove mismatches for testing
seeds <- semi_join(seeds, trust19Lookup, by = c("code" = "trustId")) %>%
    mutate(trust = as.numeric(as.factor(code)))

## write lookup table mapped to three weeks before
mutate(seeds, date = date - 21) %>%
    select(date, trust, deaths) %>%
    mutate(deaths = 100 * deaths) %>%
    write_csv("../../inputs/time_seeds.csv", col_names = FALSE)

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
seeds <- select(seeds, code, trust) %>%
    distinct() %>%
    inner_join(Ward19Lookup, by = c("code" = "trustId")) %>%
    select(FID, trust, prop) %>%
    arrange(trust, FID)
stopifnot(
    group_by(seeds, trust) %>%
    summarise(prop = sum(prop)) %>%
    pluck("prop") %>%
    {all(. == 1)}
)

## write seeds file
write_csv(seeds, "../../inputs/ward_seeds.csv", col_names = FALSE)
    