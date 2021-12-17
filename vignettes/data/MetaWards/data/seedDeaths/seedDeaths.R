## load libraries
library(tidyverse)
library(patchwork)
library(sf)
library(lubridate)
library(magrittr)

## extract command line arguments
args <- commandArgs(TRUE)
if(length(args) > 0) {
    stopifnot(length(args) == 2)
    startdate <- args[1]
    enddate <- args[2]
} else {
    stop("No arguments")
}

# ## set dates
# startdate <- "2020-02-01"
# enddate <- "2020-03-11"

## set path to MetaWardsData
path <- paste0("../../../../../../MetaWardsData/model_data/2011to2019Data/")

## read in data, originally from here: 
## https://api.coronavirus.data.gov.uk/v2/data?areaType=ltla&metric=cumDeathsByDeathDate&format=csv
deaths <- read_csv("ltla_2021-05-18.csv", guess_max = 30000) %>%
    filter(date <= enddate)

## check unique date / LAD entries
group_by(deaths, date, areaCode) %>%
    count() %>%
    pluck("n") %>%
    summary()

## load LAD to region lookup table
regionLookup <- read_csv("Local_Authority_District_to_Region_(April_2019)_Lookup_in_England.csv")

## check everything matches
temp <- anti_join(deaths, regionLookup, by = c("areaCode" = "LAD19CD"))
temp

## examine cases by date
p1 <- deaths %>%
    group_by(date) %>%
    summarise(deaths = sum(cumDeathsByDeathDate), .groups = "drop") %>%
    arrange(date) %>%
    mutate(deaths = c(deaths[1], diff(deaths))) %>%
    ggplot(aes(x = date, y = deaths)) +
    geom_line() + xlab("Date") + ylab("Deaths")
p2 <- deaths %>%
    group_by(date) %>%
    summarise(deaths = sum(cumDeathsByDeathDate), .groups = "drop") %>%
    ggplot(aes(x = date, y = deaths)) +
    geom_line() + xlab("Date") + ylab("Cumulative Deaths")
p3 <- deaths %>%
    inner_join(regionLookup, by = c("areaCode" = "LAD19CD")) %>%
    group_by(date, RGN19NM) %>%
    summarise(deaths = sum(cumDeathsByDeathDate), .groups = "drop") %>%
    arrange(date) %>%
    group_by(RGN19NM) %>%
    mutate(deaths = c(deaths[1], diff(deaths))) %>%
    ggplot(aes(x = date, y = deaths, colour = RGN19NM)) +
    geom_line() + xlab("Date") + ylab("Deaths") +
    labs(colour = "Region")
p4 <- deaths %>%
    inner_join(regionLookup, by = c("areaCode" = "LAD19CD")) %>%
    group_by(date, RGN19NM) %>%
    summarise(deaths = sum(cumDeathsByDeathDate), .groups = "drop") %>%
    ggplot(aes(x = date, y = deaths, colour = RGN19NM)) +
    geom_line() + xlab("Date") + ylab("Cumulative Deaths") +
    labs(colour = "Region")
p <- (p1 + p2) / (p3 + p4) + plot_layout(guides = "collect")
ggsave("deathCounts.pdf", p, width = 10, height = 8)

## proportions by region
p <- deaths %>%
    inner_join(regionLookup, by = c("areaCode" = "LAD19CD")) %>%
    group_by(date, RGN19NM) %>%
    summarise(deaths = sum(cumDeathsByDeathDate), .groups = "drop_last") %>%
    mutate(prop = deaths / sum(deaths)) %>%
    ggplot(aes(x = date, y = prop, fill = RGN19NM)) +
    geom_bar(stat = "identity") + xlab("Date") + 
    ylab("Proportion cumulative deaths") +
    labs(fill = "Region")
ggsave("deathsProps.pdf", p)

## now read in ward to LAD lookup
Ward19Lookup <- read_csv(paste0(path, "Ward19_Lookup.csv"))

## extract cumulative deaths before end date
seeds <- deaths %>%
    rename(cumDeaths = cumDeathsByDeathDate) %>%
    select(areaCode, date, cumDeaths) %>%
    filter(date <= enddate) %>%
    filter(deaths > 0) %>%
    complete(date = seq(ymd(startdate), ymd(enddate), by = 1), areaCode = unique(Ward19Lookup$LAD19CD)) %>%
    arrange(areaCode, date) %>%
    group_by(areaCode) %>%
    mutate(cumDeaths = ifelse(date == ymd(startdate) & is.na(cumDeaths), 0, cumDeaths)) %>%
    fill(cumDeaths) %>%
    ungroup()

## check
group_by(seeds, areaCode) %>%
    slice(1:3)

# ## quick plot
# ggplot(seeds) +
#     geom_line(aes(x = date, y = cumDeaths, group = areaCode))
    
## check LTLA names match
temp <- anti_join(seeds, Ward19Lookup, by = c("areaCode" = "LAD19CD"))
temp

## match to wards
seeds <- mutate(seeds, lad = as.numeric(as.factor(areaCode))) %>%
    group_by(areaCode, lad) %>%
    nest()

## load workers and players
workers <- read_delim(paste0(path, "WorkSize19.dat"), col_names = FALSE, delim = " ")
players <- read_delim(paste0(path, "PlaySize19.dat"), col_names = FALSE, delim = " ")
popsize <- full_join(workers, players, by = "X1") %>%
    replace_na(list(X2.x = 0, X2.y = 0)) %>%
    mutate(popsize = X2.x + X2.y) %>%
    select(-X2.x, -X2.y)

## join to Ward19 lookup
Ward19Lookup <- Ward19Lookup %>%
    inner_join(popsize, by = c("FID" = "X1")) %>%
    group_by(LAD19CD) %>%
    mutate(prop = popsize / sum(popsize)) %>%
    ungroup()

## get LAD population sizes for seeding
seeds <- group_by(Ward19Lookup, LAD19CD) %>%
    summarise(popsize = sum(popsize)) %>%
    inner_join(seeds, by = c("LAD19CD" = "areaCode")) %>%
    mutate(deaths = map_dbl(data, ~max(.$cumDeaths)))

## write seeding file
saveRDS(seeds, "../../inputs/seedsInfo.rds")

## spatial plot of seeds
lad19 <- st_read("../ward11toWard19Mapping/Local_Authority_Districts_(December_2019)_Boundaries_UK_BFC/Local_Authority_Districts_(December_2019)_Boundaries_UK_BFC.shp")

## check LTLA names match
temp <- anti_join(seeds, lad19, by = c("LAD19CD" = "lad19cd"))
temp

## join seeds to shapefile
lad19 <- left_join(lad19, mutate(seeds, deaths = map_dbl(data, ~max(.$cumDeaths))), by = c("lad19cd" = "LAD19CD"))

## plot seeds spatially
p <- mutate(lad19, deaths = ifelse(deaths == 0, NA, deaths)) %>%
    mutate(deaths = as.character(deaths)) %>%
    ggplot() + 
    geom_sf(aes(fill = deaths))
ggsave("spatialSeeds.pdf", p, height = 7, width = 3.5)

## plot seeds spatially
p <- filter(lad19, deaths > 0) %>%
    mutate(deaths = as.character(deaths)) %>%
    ggplot() + 
    geom_sf(aes(fill = deaths))
ggsave("spatialSeedsZoom.pdf", p, height = 7, width = 7)

## plot sizes by death status
p <- ggplot(seeds) +
    geom_histogram(aes(x = popsize)) +
    facet_wrap(~deaths)
ggsave("popsize.pdf", p)

## subset by LADs with initial infections
seeds <- select(seeds, LAD19CD, lad) %>%
    distinct() %>%
    inner_join(Ward19Lookup, by = "LAD19CD") %>%
    select(FID, lad, prop) %>%
    arrange(lad, FID)
stopifnot(
    group_by(seeds, lad) %>%
    summarise(prop = sum(prop), .groups = "drop") %>%
    pluck("prop") %>%
    {all(. == 1)}
)

## write seeds file
write_csv(seeds, "../../inputs/ward_seeds.csv", col_names = FALSE)
